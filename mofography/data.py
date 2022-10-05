# """
# Gathering, analyzing and maniupulating large data for mofography
# ██████╗  █████╗ ████████╗ █████╗ 
# ██╔══██╗██╔══██╗╚══██╔══╝██╔══██╗
# ██║  ██║███████║   ██║   ███████║
# ██║  ██║██╔══██║   ██║   ██╔══██║
# ██████╔╝██║  ██║   ██║   ██║  ██║
# ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝
                                

# """

# TODO figure out how to use these routines to pull data from a python file
# TODO without using boiler plate code as much as possible

async def fetch(url, session):
        async with session.get(url) as response:
            # print(await response.read())
            return await response.json()

async def bound_fetch(sem, url, session):
        # Getter function with semaphoe.
        async with sem:
            return await fetch(url, session)


async def get_mofdb_data(session, url,params, headers, npages, limit=1000):
    """Compile a asyncio Future containing the results of a search by making requests to MOFDb in parallel asynchronously

    :param session: object of the ClientSession
    :type session: aiohttp.ClientSession
    :param url: URL of the website. Designed for MOFDb, so may not work with other URLs, defaults to 'https://mof.tech.northwestern.edu/mofs/#.json'
    :type url: str, optional
    :param headers: dictionary of search parameters, defaults to {'loading': 'cm3(STP)/cm3', 'pressure': 'bar'}
    :type headers: dict, optional
    :param params: dictionary of the headers to be output in the result, defaults to {'vf_min': 0.75, 'vf_max': 1, 'database':'COREMOF 2019','page':1}
    :type params: dict, optional
    :param npages: number of the result pages to be output, defaults to 1
    :type npages: int, optional
    :param limit: semaphore limit of number of jobs per block, no need to change this, defaults to 1000
    :type limit: int, optional
    :return: an instance of asyncio.Future that contains all the search results
    :rtype: asyncio.Future
    """        
    from furl import furl
    import asyncio
    tasks = []
    # create instance of Semaphore
    sem = asyncio.Semaphore(limit)
    for i in range(1,npages+1):

        url_string = furl(url)
        params['page']=i
        url_string.set(params)

        # pass Semaphore and session to every GET request
        task = asyncio.ensure_future(bound_fetch(sem, url_string.tostr(), session))
        tasks.append(task)

    responses = asyncio.gather(*tasks)
    return await responses

def get_number_of_pages(url, headers, params):

    """
    Get the number of pages of the result, from the url, for a given set of headers and parameters

    :param url: url for MOFDb for example
    :type url: string
    :param headers: dictionary of headers
    :type headers: dict
    :param params: dictionary of parameters
    :type params: dict
    :return: number of pages in the result
    :rtype: int
    """    

    import requests
    resp = requests.get(url=url, headers=headers, params=params)
    resp = resp.json()
    return resp['pages']

def get_data_from_mofdb(url='https://mof.tech.northwestern.edu/mofs/#.json', params={'vf_min': 0.1, 'vf_max': 1.0, 'database':'IZA','page':1} , headers={'loading': 'cm3(STP)/cm3', 'pressure': 'bar'},npages=None, limit=1000):
    """
    This is a wrapper for the MOFDB API, which is explained here:https://mof.tech.northwestern.edu/api 
    It takes in a url, parameters, headers, number of pages and a path to the mofography repository and
    returns a pandas dataframe of the data. 
    
    The function is asynchronous and uses the aiohttp library to make multiple requests to the MOFDB
    API. 
    
    
    
    :param url: The url of the MOF database, defaults to https://mof.tech.northwestern.edu/mofs/#.json
    (optional)
    :param params: a dictionary of parameters to be passed to the url. The default parameters are: {'vf_min': 0.1, 'vf_max': 1.0, 'database':'IZA','page':1}
    :param headers: This is the headers that you want to pass to the MOFDB API. The default is to get
    the loading and pressure data
    :param npages: number of pages to be scraped. If None, then the number of pages is determined by the
    function to return all the pages
    :param limit: maximum number of pages allowed in the query (no need to change this one unless overflow )
    :return: A dataframe with the results of the query.
    """

    from aiohttp import ClientSession
    import asyncio
    import numpy as np
    import nest_asyncio
    nest_asyncio.apply()
    import pandas as pd 


    if npages==None:
        npages=get_number_of_pages(url, headers, params)
    
    # This shouldn't be called in a function (similar to initiating a Dask cluster)
    loop = asyncio.get_event_loop()
    # with ClientSession(headers= headers) as session:
    session = ClientSession(headers=headers)
    fut = asyncio.ensure_future(get_mofdb_data(session, url=url, params=params, headers=headers, npages=npages, limit=limit))
    loop.run_until_complete(fut)
    session.close()
    data_df = pd.DataFrame(list(np.hstack(([f['results'] for f in fut.result()]))))

    return data_df   

def atoms_from_mofdb_df(mofdb_df):

    """
    This function takes a MOFDB dataframe from get_data_from_mofdb routine, which has a column called 'cif' that contains CIF data as strings, and
    returns a list of ASE atoms objects
    
    :param mofdb_df: the dataframe containing the MOFDB data
    :return: A list of atoms objects.
    """
    import numpy as np
    def read_cif_string(cif):
        """
        It takes a string of CIF data and returns an ASE atoms object.
        
        :param cif: The CIF file as a string
        :return: the ASE atoms object.
        """
        from ase.io import read 
        import io
        try:
            output = read(io.BytesIO(bytes(cif, encoding='utf-8')), format='cif')
        except: 
            output = 'NaN'
        return output

    # from rich.progress import track
    atoms_list = [read_cif_string(x) for x in mofdb_df['cif']]
    return np.array(atoms_list)

def read_raspa_pdb(path_to_file):
        """
            created by Arun Gopalan, Snurr Research Group.                   _                                 ___  ___  ___ 
            _ __ ___  __ _  __| |  _ __ __ _ ___ _ __   __ _    / _ \/   \/ __\
            | '__/ _ \/ _` |/ _` | | '__/ _` / __| '_ \ / _` |  / /_)/ /\ /__\//
            | | |  __/ (_| | (_| | | | | (_| \__ \ |_) | (_| | / ___/ /_// \/  \
            |_|  \___|\__,_|\__,_| |_|  \__,_|___/ .__/ \__,_| \/  /___,'\_____/
                                                |_|                            
            Reads the output PDB movie file from RASPA into separate configurations, including the chemical symbols and
            cell dimensions. 
            :type path_to_file: str
            :param path_to_file: path to the RASPA PDB file
        
            :raises:
            Not sure.
        
            :rtype: a python dictionary with 'cells', 'symbols' and 'coord' which stand for cell dimensions, printed symbols and coordinates
            """

        import numpy as np

        f = open(path_to_file).readlines()
        start = np.where(["MODEL" in line for line in f])[0] + 2  # * Start of config
        ends = np.where(["ENDMDL" in line for line in f])[0]  # End for config
        cryst = np.where(["CRYST" in line for line in f])[0]  # box shape for the config

        data = [f[start[i]:ends[i]] for i in range(len(start))]

        coord = np.array([[np.array(line.split()[4:7]).astype(float) for line in d] for d in data])
        cell_dims = np.array([np.array(line.split())[1:].astype(float) for line in f if "CRYST" in line])
        symbols = np.array([[np.array(line.split()[2]) for line in d] for d in data])

        output = {}
        output['cells'] = cell_dims
        output['coords'] = coord
        output['symbols'] = symbols

        return output