"""
Gathering, analyzing and maniupulating large data for mofography
██████╗  █████╗ ████████╗ █████╗ 
██╔══██╗██╔══██╗╚══██╔══╝██╔══██╗
██║  ██║███████║   ██║   ███████║
██║  ██║██╔══██║   ██║   ██╔══██║
██████╔╝██║  ██║   ██║   ██║  ██║
╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝
                                

"""
async def fetch(url, session):
        async with session.get(url) as response:
            # print(await response.read())
            return await response.json()

async def bound_fetch(sem, url, session):
        # Getter function with semaphoe.
        async with sem:
            return await fetch(url, session)


async def get_mofdb_data(session, url='https://mof.tech.northwestern.edu/mofs/#.json', headers={'loading': 'cm3(STP)/cm3', 'pressure': 'bar'}, params={'vf_min': 0.75, 'vf_max': 1, 'database':'COREMOF 2019','page':1}, npages=1, limit=1000):
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