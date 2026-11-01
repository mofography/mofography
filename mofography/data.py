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
    """Fetch JSON from a URL with aiohttp (single-page MOFdb helper)."""
    async with session.get(url) as response:
        # print(await response.read())
        return await response.json()


async def bound_fetch(sem, url, session):
    """Fetch JSON from a URL while holding an asyncio semaphore (parallel MOFdb page downloads)."""
    async with sem:
        return await fetch(url, session)


async def get_mofdb_data(session, url, params, headers, npages, limit=1000):
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
    from rich.progress import track

    # create instance of Semaphore
    sem = asyncio.Semaphore(limit)
    for i in track(range(1, npages + 1)):
        url_string = furl(url)
        params["page"] = i
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
    return resp["pages"]


def create_db_from_pickled_graphs(db_path, pickle_path):
    """
    Create a SQLite database and store gzip-compressed pickled graphs from a pickle file.

    :param db_path: The `db_path` parameter is the path to the SQLite database where you want to store
    the graphs. This is the location where the graphs will be saved after being pickled, compressed, and
    inserted into the database
    :param pickle_path: The `pickle_path` parameter in the `create_db_from_pickled_graphs` function is
    the file path to the pickled collection of graphs that you want to store in the SQLite database.
    This file contains the serialized graphs that need to be unpickled, compressed, and stored in the
    database
    """

    import sqlite3
    import pickle
    import gzip

    # Connect to the SQLite database and create a table to store graphs
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS graphs
                 (id INTEGER PRIMARY KEY, data BLOB)""")

    # Unpickle your collection of graphs
    with open(pickle_path, "rb") as f:
        graphs = pickle.load(f)

    # For each graph, pickle it, compress it, and store it in the SQLite database
    for i, G in enumerate(graphs):
        pickled_graph = pickle.dumps(G)
        compressed_graph = gzip.compress(pickled_graph)
        c.execute("INSERT INTO graphs (data) VALUES (?)", (compressed_graph,))
    conn.commit()


def get_graph(db_path, id):
    """
    Load one gzip-pickled graph from a SQLite graphs table by row id.

    :param db_path: The `db_path` parameter is the file path to the SQLite database where the graph data
    is stored. This function connects to the SQLite database at the specified path to retrieve the graph
    data
    :param id: The `id` parameter in the `get_graph` function is used to specify which graph to retrieve
    from the database. It is likely a unique identifier that helps identify a specific graph stored in
    the database table. When calling the `get_graph` function, you would provide the `id` of the
    :return: The function `get_graph` returns a graph object that is retrieved from a SQLite database.
    The graph is decompressed using gzip, unpickled, and then returned as the final output.
    """

    import sqlite3
    import pickle
    import gzip

    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Retrieve the graph from the database, decompress it, unpickle it, and return it
    c.execute("SELECT data FROM graphs WHERE id=?", (id,))
    compressed_graph = c.fetchone()[0]
    pickled_graph = gzip.decompress(compressed_graph)
    return pickle.loads(pickled_graph)


def add_graph_to_db(db_path, graph):
    """
    Insert one gzip-pickled NetworkX (or similar) graph into a SQLite graphs table.

    :param db_path: The `db_path` parameter is the file path to the SQLite database where you want to
    store the graph data. This function connects to the SQLite database at the specified path and
    inserts the compressed graph data into the database
    :param graph: The `graph` parameter in the `add_graph_to_db` function is the graph data structure
    that you want to store in the SQLite database. This graph can be represented as a Python object
    (e.g., dictionary, list of edges, etc.) that you want to save for later retrieval or analysis
    """

    import sqlite3
    import pickle
    import gzip

    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Pickle the graph, compress it, and store it in the SQLite database
    pickled_graph = pickle.dumps(graph)
    compressed_graph = gzip.compress(pickled_graph)
    c.execute("INSERT INTO graphs (data) VALUES (?)", (compressed_graph,))
    conn.commit()


def get_num_graphs(db_path):
    """
    Return the number of rows in the SQLite ``graphs`` table at ``db_path``.

    :param db_path: The `db_path` parameter in the `get_num_graphs` function is the path to the SQLite
    database file that contains the table of graphs. This function connects to the SQLite database at
    the specified path, executes a SQL query to count the number of rows in the `graphs` table, and
    :return: The function `get_num_graphs` returns the number of rows in the "graphs" table from the
    SQLite database specified by the `db_path`.
    """

    import sqlite3

    # Connect to the SQLite database
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Execute a SQL query that counts the number of rows in the table
    c.execute("SELECT COUNT(*) FROM graphs")
    num_graphs = c.fetchone()[0]

    return num_graphs


async def get_data_from_mofdb(
    url="https://mof.tech.northwestern.edu/mofs.json",
    params={"vf_min": 0.1, "vf_max": 1.0, "database": "IZA", "page": 1},
    headers={"loading": "cm3(STP)/cm3", "pressure": "bar"},
    npages=None,
    limit=1000,
):
    """Download MOFdb search pages in parallel and return a combined results DataFrame."""
    import aiohttp
    import pandas as pd
    import numpy as np

    # copy so we don't mutate the caller's dict; page is set for pagination
    params = dict(params)
    params["page"] = 1

    if npages is None:
        npages = get_number_of_pages(url, headers, params)

    # Use a very generous timeout since some mofdb pages can be slow
    timeout = aiohttp.ClientTimeout(total=None, sock_read=None)
    async with aiohttp.ClientSession(headers=headers, timeout=timeout) as session:
        results = await get_mofdb_data(
            session, url=url, params=params, headers=headers, npages=npages, limit=limit
        )
    data_df = pd.DataFrame(list(np.hstack(([f["results"] for f in results]))))
    return data_df


# Default columns to stratify on for representative sampling (MOFdb-style df)
MOFDB_STRAT_COLS = ["surface_area_m2g", "void_fraction", "lcd", "pld"]


def sample_representative_mofdb(
    df,
    n,
    strat_cols=None,
    n_bins=4,
    random_state=0,
):
    """
    Sample a representative set of N materials from a MOFdb-style dataframe using
    proportional stratified sampling so the sample's marginal distributions match
    the population over the given columns (e.g. surface area, void fraction, LCD, PLD).

    :param df: DataFrame with MOFdb-style columns (e.g. surface_area_m2g, void_fraction, lcd, pld).
    :param n: Target number of materials to sample.
    :param strat_cols: List of column names to stratify on. Defaults to MOFDB_STRAT_COLS.
    :param n_bins: Number of quantile-based bins per column for stratification.
    :param random_state: Random seed for reproducibility.
    :return: DataFrame of sampled rows (materials).
    """
    import numpy as np
    import pandas as pd

    strat_cols = strat_cols or MOFDB_STRAT_COLS
    missing = [c for c in strat_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"DataFrame missing columns: {missing}. Available: {list(df.columns)}"
        )
    df = df.dropna(subset=strat_cols)
    N = len(df)
    if n >= N:
        return df.sample(n=min(n, N), random_state=random_state)

    stratum_parts = []
    for c in strat_cols:
        labels, _ = pd.qcut(
            df[c], q=n_bins, labels=False, retbins=True, duplicates="drop"
        )
        stratum_parts.append(labels.astype(str))
    df = df.copy()
    df["_stratum"] = stratum_parts[0]
    for s in stratum_parts[1:]:
        df["_stratum"] = df["_stratum"] + "_" + s

    g = df.groupby("_stratum", group_keys=False)
    counts = g.size()
    exact = n * counts / counts.sum()
    floor = np.floor(exact).astype(int)
    remainder = exact - floor
    to_take = floor.copy()
    shortfall = n - int(to_take.sum())
    if shortfall > 0:
        add_idx = remainder.nlargest(min(shortfall, len(remainder))).index
        to_take.loc[add_idx] += 1

    sampled = []
    for stratum_id in to_take.index:
        size = int(to_take.loc[stratum_id])
        if size <= 0:
            continue
        block = g.get_group(stratum_id)
        k = min(size, len(block))
        sampled.append(block.sample(n=k, random_state=random_state))
    out = pd.concat(sampled, axis=0)
    return out.drop(columns=["_stratum"], errors="ignore")


def plot_representative_sample_subplots(
    full,
    sample,
    strat_cols=None,
    nbins_hist=40,
    n_sample=None,
    title=None,
):
    """
    Build a 2x2 Plotly figure comparing full vs sample histograms for each
    stratification column (surface area, void fraction, LCD, PLD). Use to verify
    that a representative sample matches the full distribution.

    :param full: Full population DataFrame (MOFdb-style columns).
    :param sample: Sampled DataFrame from sample_representative_mofdb.
    :param strat_cols: Columns to plot; defaults to MOFDB_STRAT_COLS.
    :param nbins_hist: Number of histogram bins per subplot.
    :param n_sample: Number to show in title (default: len(sample)).
    :param title: Main figure title. If None, uses "Representative subset (N={n_sample})".
    :return: Plotly Figure.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    strat_cols = strat_cols or MOFDB_STRAT_COLS
    n_sample = n_sample if n_sample is not None else len(sample)
    labels = {
        "surface_area_m2g": "Surface area (m²/g)",
        "void_fraction": "Void fraction",
        "lcd": "LCD (Å)",
        "pld": "PLD (Å)",
    }
    if title is None:
        title = f"Representative subset (N={n_sample})"

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("", "", "", ""),
        vertical_spacing=0.12,
        horizontal_spacing=0.18,
    )
    for i, col in enumerate(strat_cols):
        if col not in full.columns or col not in sample.columns:
            continue
        row, c = i // 2 + 1, i % 2 + 1
        x_full = full[col].dropna()
        x_samp = sample[col].dropna()
        low = min(x_full.min(), x_samp.min())
        high = max(x_full.max(), x_samp.max())
        if low == high:
            low, high = low - 0.5, high + 0.5
        xbins = dict(
            start=float(low),
            end=float(high),
            size=float((high - low) / nbins_hist),
        )
        show_leg = i == 0
        fig.add_trace(
            go.Histogram(
                x=x_full,
                xbins=xbins,
                name="Full (>0.2 Vf)",
                marker=dict(
                    color="rgba(100,100,100,0.5)", line=dict(width=0.5)
                ),
                legendgroup="full",
                showlegend=show_leg,
            ),
            row=row,
            col=c,
        )
        fig.add_trace(
            go.Histogram(
                x=x_samp,
                xbins=xbins,
                name="Sample",
                marker=dict(
                    color="rgba(30,120,180,0.7)", line=dict(width=0.5)
                ),
                legendgroup="sample",
                showlegend=show_leg,
            ),
            row=row,
            col=c,
        )
        fig.update_xaxes(title_text=labels.get(col, col), row=row, col=c)
        fig.update_yaxes(title_text="Count", row=row, col=c)

    font_arial_16 = dict(family="Arial", size=16)
    fig.update_layout(
        template="plotly_white",
        title=dict(text=title, x=0.5, font=dict(family="Arial", size=16)),
        barmode="overlay",
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5,
            font=font_arial_16,
        ),
        margin=dict(t=80, b=55, l=50, r=40),
        height=500,
        width=700,
    )
    fig.update_xaxes(title_font=font_arial_16)
    fig.update_yaxes(title_font=font_arial_16)
    axes_style = dict(
        showline=True,
        linecolor="black",
        mirror=True,
        zeroline=True,
        zerolinecolor="black",
        rangemode="tozero",
    )
    fig.update_xaxes(**axes_style)
    fig.update_yaxes(**axes_style)
    return fig


# def get_data_from_mofdb(
#     url="https://mof.tech.northwestern.edu/mofs/#.json",
#     params={"vf_min": 0.1, "vf_max": 1.0, "database": "IZA", "page": 1},
#     headers={"loading": "cm3(STP)/cm3", "pressure": "bar"},
#     npages=None,
#     limit=1000,
# ):
#     """
#     This is a wrapper for the MOFDB API, which is explained here:https://mof.tech.northwestern.edu/api
#     It takes in a url, parameters, headers, number of pages and a path to the mofography repository and
#     returns a pandas dataframe of the data.
#
#     The function is asynchronous and uses the aiohttp library to make multiple requests to the MOFDB
#     API.
#
#
#
#     :param url: The url of the MOF database, defaults to https://mof.tech.northwestern.edu/mofs/#.json
#     (optional)
#     :param params: a dictionary of parameters to be passed to the url. The default parameters are: {'vf_min': 0.1, 'vf_max': 1.0, 'database':'IZA','page':1}
#     :param headers: This is the headers that you want to pass to the MOFDB API. The default is to get
#     the loading and pressure data
#     :param npages: number of pages to be scraped. If None, then the number of pages is determined by the
#     function to return all the pages
#     :param limit: maximum number of pages allowed in the query (no need to change this one unless overflow )
#     :return: A dataframe with the results of the query.
#     """
#
#     from aiohttp import ClientSession
#     import asyncio
#     import numpy as np
#     import nest_asyncio
#
#     nest_asyncio.apply()
#     import pandas as pd
#
#     # just put page:1 in the params it doesn't matter, only the npages parameter is used
#     params["page"] = 1
#
#     if npages == None:
#         npages = get_number_of_pages(url, headers, params)
#
#     # This shouldn't be called in a function (similar to initiating a Dask cluster)
#     loop = asyncio.get_event_loop()
#     # with ClientSession(headers= headers) as session:
#     session = ClientSession(headers=headers)
#     fut = asyncio.ensure_future(
#         get_mofdb_data(
#             session, url=url, params=params, headers=headers, npages=npages, limit=limit
#         )
#     )
#     loop.run_until_complete(fut)
#     session.close()
#     data_df = pd.DataFrame(list(np.hstack(([f["results"] for f in fut.result()]))))
#
#     return data_df


def read_cif_string(cif):
    """
    Parse a CIF string into an ASE Atoms object, or the string ``NaN`` on failure.

    :param cif: The CIF file as a string
    :return: the ASE atoms object.
    """
    from ase.io import read
    import io

    try:
        output = read(io.BytesIO(bytes(cif, encoding="utf-8")), format="cif")
    except:
        output = "NaN"
    return output


def atoms_from_mofdb_df(mofdb_df):
    """
    Parse each MOFdb DataFrame ``cif`` column to ASE Atoms and attach ``info['name']`` from ``name``.

    :param mofdb_df: the dataframe containing the MOFDB data
    :return: A list of atoms objects.
    """
    import numpy as np
    from tqdm import tqdm

    # import dask.bag as db
    # from dask.diagnostics import ProgressBar
    # with ProgressBar():
    #     atoms_list = db.from_sequence(mofdb_df['cif'], npartitions=20).map(read_cif_string).compute()
    # from rich.progress import track
    atoms_list = [read_cif_string(x) for x in tqdm(mofdb_df["cif"])]
    for i, a in enumerate(atoms_list):
        if a != "NaN":
            a.info["name"] = mofdb_df["name"][i]

    return np.array(atoms_list, dtype="object")


def get_atoms_from_mofdb_df_parallel(
    mofdb_df,
    scheduler: str = "threads",
):
    """
    Parallel variant of ``atoms_from_mofdb_df`` for large MOFdb-style DataFrames.

    Uses Dask delayed to parse CIF strings into ASE Atoms objects in parallel.

    - Input df is expected to have columns:
        - ``'cif'`` : CIF content as strings
        - ``'name'``: material name (used to set ``atoms.info['name']``)
    - Each failed or malformed CIF is turned into the string ``'NaN'``,
      matching the behavior of ``read_cif_string`` and ``atoms_from_mofdb_df``.

    Parameters
    ----------
    mofdb_df : pandas.DataFrame
        MOFdb-style dataframe from ``get_data_from_mofdb`` (or similar) with
        at least ``'cif'`` and ``'name'`` columns.
    scheduler : str, optional
        Dask scheduler to use for ``dask.compute`` (default: ``'threads'``).

    Returns
    -------
    numpy.ndarray
        1D object array of ASE Atoms objects or ``'NaN'`` entries, aligned
        with the rows of ``mofdb_df``.
    """
    import numpy as np
    from dask import delayed, compute

    if "cif" not in mofdb_df.columns:
        raise ValueError("DataFrame missing required column 'cif'")
    if "name" not in mofdb_df.columns:
        raise ValueError("DataFrame missing required column 'name'")

    def _safe_read_atoms(cif_str, name):
        atoms = read_cif_string(cif_str)
        if atoms != "NaN":
            atoms.info["name"] = name
        return atoms

    names = mofdb_df["name"].tolist()
    cif_strings = mofdb_df["cif"].tolist()

    delayed_atoms = [
        delayed(_safe_read_atoms)(cif, name)
        for cif, name in zip(cif_strings, names)
    ]

    atoms_list = list(compute(*delayed_atoms, scheduler=scheduler))
    return np.array(atoms_list, dtype="object")


def write_cifs_from_mofdb_df(
    mofdb_df,
    folder,
    name_column="name",
    cif_column="cif",
    skip_failed=True,
):
    """
    Write CIF strings from a MOFdb-style DataFrame to .cif files in a folder.
    Creates the folder if it does not exist. Uses the name column for filenames
    (sanitized); duplicate names get _1, _2, ... appended.

    :param mofdb_df: DataFrame with CIF data (e.g. from get_data_from_mofdb or
        sample_representative_mofdb). Must have columns for names and CIF strings.
    :param folder: Path or str to the output directory (created if missing).
    :param name_column: Column used for filenames (default "name").
    :param cif_column: Column containing CIF content as strings (default "cif").
    :param skip_failed: If True, skip rows with missing or empty CIF (default True).
    :return: List of Paths of written .cif files.
    """
    import re
    from pathlib import Path

    folder = Path(folder).expanduser().resolve()
    folder.mkdir(parents=True, exist_ok=True)

    if name_column not in mofdb_df.columns:
        raise ValueError(f"DataFrame missing column: {name_column}")
    if cif_column not in mofdb_df.columns:
        raise ValueError(f"DataFrame missing column: {cif_column}")

    def safe_filename(s):
        s = str(s).strip()
        s = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", s)
        return s[:200] if s else "unnamed"

    written = []
    seen = {}
    for idx, row in mofdb_df.iterrows():
        cif = row[cif_column]
        if skip_failed and (cif is None or (isinstance(cif, str) and not cif.strip())):
            continue
        if not isinstance(cif, str):
            cif = str(cif)
        base = safe_filename(row[name_column])
        stem = base
        suffix = ""
        n = 0
        while True:
            path = folder / f"{stem}{suffix}.cif"
            if path not in written and not path.exists():
                break
            n += 1
            suffix = f"_{n}"
            stem = base
        path.write_text(cif, encoding="utf-8")
        written.append(path)
    return written


def read_graspa_lammps_data(path_to_file):
    """Read a gRASPA / RASPA ``WriteLammpsdata`` LAMMPS ``.data`` snapshot.

    Parses the header (orthogonal or triclinic box), ``Masses`` type labels,
    and ``Atoms # charge`` records written by gRASPA movie output
    (``Movies/System_0/result_*.data``).

    Each atom record stores:

    - ``id``, ``mol_id``, ``type_id``, ``charge``, ``position`` (3,)
    - ``mass_label`` from the ``Masses`` section (e.g. ``CH3_xyl``, ``C``)
    - ``comment`` trailing text after ``#``
    - ``molecule_name`` when present in the comment (e.g. ``o-xylene``)
    - ``symbol`` element tag for ASE (framework from ``.cif`` comment token;
      guests from mass label with united-atom xylene types mapped to C/H)

    Returns
    -------
    dict
        ``n_atoms``, ``cell`` (3×3), ``box`` (xlo/xhi/...), ``type_labels``,
        ``atoms`` (list of atom dicts), ``source_path``.
    """
    import re
    from pathlib import Path

    import numpy as np

    path = Path(path_to_file).expanduser().resolve()
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    n_atoms = None
    box = {}
    type_labels: dict[int, str] = {}
    atoms: list[dict] = []

    section = None
    mass_re = re.compile(r"^\s*(\d+)\s+([0-9.eE+-]+)\s*(?:#\s*(.+))?\s*$")

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        if line.endswith(" atoms"):
            n_atoms = int(line.split()[0])
            section = None
            continue

        if line.endswith(" atom types"):
            section = None
            continue

        if line.endswith(" xlo xhi"):
            parts = line.split()
            box["xlo"], box["xhi"] = float(parts[0]), float(parts[1])
            continue
        if line.endswith(" ylo yhi"):
            parts = line.split()
            box["ylo"], box["yhi"] = float(parts[0]), float(parts[1])
            continue
        if line.endswith(" zlo zhi"):
            parts = line.split()
            box["zlo"], box["zhi"] = float(parts[0]), float(parts[1])
            continue
        if line.endswith(" xy xz yz"):
            parts = line.split()
            box["xy"], box["xz"], box["yz"] = (
                float(parts[0]),
                float(parts[1]),
                float(parts[2]),
            )
            continue

        if line == "Masses":
            section = "masses"
            continue
        if line.startswith("Pair ") or line.startswith("Bond ") or line.startswith("Angle "):
            section = None
            continue
        if line == "Atoms":
            section = "atoms"
            continue

        if section == "masses":
            m = mass_re.match(line)
            if m:
                tid = int(m.group(1))
                label = (m.group(3) or "").strip()
                type_labels[tid] = label
            continue

        if section == "atoms":
            comment = ""
            body = line
            if "#" in line:
                body, comment = line.split("#", 1)
                comment = comment.strip()
            toks = body.split()
            if len(toks) < 7:
                continue
            aid = int(toks[0])
            mol_id = int(toks[1])
            type_id = int(toks[2])
            charge = float(toks[3])
            pos = np.array([float(toks[4]), float(toks[5]), float(toks[6])])
            mass_label = type_labels.get(type_id, f"type_{type_id}")
            molecule_name = None
            symbol = None
            if comment:
                parts = comment.split()
                if parts and parts[0] in ("o-xylene", "m-xylene", "p-xylene"):
                    molecule_name = parts[0]
                if ".cif" in comment and len(parts) >= 2:
                    symbol = parts[-1]
            if symbol is None:
                symbol = _graspa_mass_label_to_symbol(mass_label)
            atoms.append(
                {
                    "id": aid,
                    "mol_id": mol_id,
                    "type_id": type_id,
                    "charge": charge,
                    "position": pos,
                    "mass_label": mass_label,
                    "comment": comment,
                    "molecule_name": molecule_name,
                    "symbol": symbol,
                }
            )

    if n_atoms is None:
        raise ValueError(f"Could not find atom count in {path}")
    if len(atoms) != n_atoms:
        raise ValueError(
            f"Expected {n_atoms} atoms in {path}, parsed {len(atoms)}"
        )
    required = ("xlo", "xhi", "ylo", "yhi", "zlo", "zhi", "xy", "xz", "yz")
    if any(k not in box for k in required):
        missing = [k for k in required if k not in box]
        raise ValueError(f"Missing box fields {missing} in {path}")

    lx = box["xhi"] - box["xlo"]
    ly = box["yhi"] - box["ylo"]
    lz = box["zhi"] - box["zlo"]
    xy, xz, yz = box["xy"], box["xz"], box["yz"]
    cell = np.array(
        [
            [lx, 0.0, 0.0],
            [xy, ly, 0.0],
            [xz, yz, lz],
        ],
        dtype=float,
    )

    return {
        "n_atoms": n_atoms,
        "cell": cell,
        "box": box,
        "type_labels": type_labels,
        "atoms": atoms,
        "source_path": str(path),
    }


def _graspa_mass_label_to_symbol(mass_label: str) -> str:
    """Map gRASPA mass-table labels to ASE element symbols."""
    label = (mass_label or "").strip()
    if not label:
        return "X"
    if label in ("H_xyl",):
        return "H"
    if label in ("C_xyl", "CH3_xyl"):
        return "C"
    # Framework labels are plain element names (Al, C, Co, ...).
    return label.split("_")[0]


def read_raspa_pdb(path_to_file, symbol_map=None):
    """
    Parse a RASPA trajectory PDB into per-frame cells, symbols, and Cartesian coordinates.
            :type path_to_file: str
            :param path_to_file: path to the RASPA PDB file
            :param symbol_map: Optional dictionary mapping custom symbols to standard element symbols.
        
            :raises:
            Not sure.
        
            :rtype: a python dictionary with 'cells', 'symbols' and 'coord' which stand for cell dimensions, printed symbols and coordinates
            """

    import numpy as np

    f = open(path_to_file).readlines()
    start = np.where(["MODEL" in line for line in f])[0] + 2  # * Start of config
    ends = np.where(["ENDMDL" in line for line in f])[0]  # End for config
    cryst = np.where(["CRYST" in line for line in f])[0]  # box shape for the config

    data = [f[start[i] : ends[i]] for i in range(len(start))]

    # Initialize lists to hold coordinates and symbols
    coords = []
    symbols = []  # Change: Initialize symbols as a list of lists
    skipped_lines = []  # List to hold skipped lines

    for d in data:
        coord = []
        config_symbols = []  # New list to hold symbols for the current configuration
        for line in d:
            parts = line.split()
            if len(parts) >= 7:  # Ensure there are enough parts to avoid index errors
                try:
                    coord.append(np.array(parts[4:7], dtype=float))
                    # Prefer the element column (last token in PDB ATOM records),
                    # fall back to the atom label when missing.
                    symbol = parts[-1] if len(parts) >= 10 else parts[2]
                    if symbol_map and symbol in symbol_map:
                        symbol = symbol_map[
                            symbol
                        ]  # Replace with mapped symbol if exists
                    config_symbols.append(
                        symbol
                    )  # Append the standardized symbol to the current config
                except ValueError:
                    skipped_lines.append(
                        line.strip()
                    )  # Add the skipped line to the list
                    continue  # Skip lines that cannot be converted to float

        coords.append(np.array(coord))
        symbols.append(
            config_symbols
        )  # Change: Append the list of symbols for the current config

    cell_dims = np.array(
        [np.array(line.split())[1:].astype(float) for line in f if "CRYST" in line]
    )

    output = {
        "cells": cell_dims,
        "coords": np.array(
            coords, dtype=object
        ),  # Use dtype=object to handle varying lengths
        "symbols": np.array(symbols, dtype=object),  # Use dtype=object for symbols
    }

    # Print warning if any lines were skipped
    if skipped_lines:
        print("Warning: The following lines were skipped due to conversion issues:")
        for skipped in skipped_lines:
            print(f"  - {skipped}")

    return output


global forcefields
forcefields = dict(
    uff={
        "Ac": [16.60, 3.10],
        "Ag": [18.11, 2.80],
        "Al": [254.09, 4.01],
        "Am": [7.04, 3.01],
        "Ar": [93.08, 3.45],
        "As": [155.47, 3.77],
        "At": [142.89, 4.23],
        "Au": [19.62, 2.93],
        "B": [90.57, 3.64],
        "Ba": [183.15, 3.30],
        "Be": [42.77, 2.45],
        "Bi": [260.63, 3.89],
        "Bk": [6.54, 2.97],
        "Br": [126.29, 3.73],
        "C": [52.83, 3.43],
        "Ca": [119.75, 3.03],
        "Cd": [114.72, 2.54],
        "Ce": [6.54, 3.17],
        "Cf": [6.54, 2.95],
        "Cl": [114.21, 3.52],
        "Cm": [6.54, 2.96],
        "Co": [7.04, 2.56],
        "Cr": [7.55, 2.69],
        "Cu": [2.52, 3.11],
        "Cs": [22.64, 4.02],
        "Dy": [3.52, 3.05],
        "Eu": [4.03, 3.11],
        "Er": [3.52, 3.02],
        "Es": [6.04, 2.94],
        "F": [25.16, 3.00],
        "Fe": [6.54, 2.59],
        "Fm": [6.04, 2.93],
        "Fr": [25.16, 4.37],
        "Ga": [208.81, 3.90],
        "Ge": [190.69, 3.81],
        "Gd": [4.53, 3.00],
        "H": [22.14, 2.57],
        "He": [10.9, 2.64],
        "Hf": [36.23, 2.80],
        "Hg": [193.71, 2.41],
        "Ho": [3.52, 3.04],
        "I": [170.57, 4.01],
        "In": [301.39, 3.98],
        "Ir": [36.73, 2.53],
        "K": [17.61, 3.40],
        "Kr": [110.69, 3.69],
        "La": [8.55, 3.14],
        "Li": [12.58, 2.18],
        "Lu": [20.63, 3.24],
        "Lr": [5.53, 2.88],
        "Md": [5.53, 2.92],
        "Mg": [55.85, 2.69],
        "Mn": [6.54, 2.64],
        "Mo": [28.18, 2.72],
        "N": [34.72, 3.26],
        "Na": [15.09, 2.66],
        "Ne": [21.13, 2.66],
        "Nb": [29.69, 2.82],
        "Nd": [5.03, 3.18],
        "No": [5.53, 2.89],
        "Ni": [7.55, 2.52],
        "Np": [9.56, 3.05],
        "O": [30.19, 3.12],
        "Os": [18.62, 2.78],
        "P": [153.46, 3.69],
        "Pa": [11.07, 3.05],
        "Pb": [333.59, 3.83],
        "Pd": [24.15, 2.58],
        "Pm": [4.53, 3.16],
        "Po": [163.52, 4.20],
        "Pr": [5.03, 3.21],
        "Pt": [40.25, 2.45],
        "Pu": [8.05, 3.05],
        "Ra": [203.27, 3.28],
        "Rb": [20.13, 3.67],
        "Re": [33.21, 2.63],
        "Rh": [26.67, 2.61],
        "Rn": [124.78, 4.25],
        "Ru": [28.18, 2.64],
        "S": [137.86, 3.59],
        "Sb": [225.91, 3.94],
        "Sc": [9.56, 2.94],
        "Se": [146.42, 3.75],
        "Si": [202.27, 3.83],
        "Sm": [4.03, 3.14],
        "Sn": [285.28, 3.91],
        "Sr": [118.24, 3.24],
        "Ta": [40.75, 2.82],
        "Tb": [3.52, 3.07],
        "Tc": [24.15, 2.67],
        "Te": [200.25, 3.98],
        "Th": [13.08, 3.03],
        "Ti": [8.55, 2.83],
        "Tl": [342.14, 3.87],
        "Tm": [3.02, 3.01],
        "U": [11.07, 3.02],
        "V": [8.05, 2.80],
        "W": [33.71, 2.73],
        "Xe": [167.04, 3.92],
        "Y": [36.23, 2.98],
        "Yb": [114.72, 2.99],
        "Zn": [62.39, 2.46],
        "Zr": [34.72, 2.78],
        "CH4_sp3": [148.0, 3.73],
        "CH3_sp3": [98.0, 3.75],
        "CH2_sp3": [46.0, 3.95],
        "CH_sp2": [50.0, 3.695],
        "CH_sp3": [10, 4.68],
        "C_sp3": [0.5, 6.4],
        "H_h2": [0.0, 0.0],
        "H_com": [36.7, 2.958],
        "N_n2": [36.0, 3.31],
        "N_com": [0.0, 0.0],
        "O_co2": [79.0, 3.05],
        "C_co2": [27.0, 2.80],
        "Ow": [78.0, 3.15],
    }
)


def mix_lorentz_berthelot(ase_atoms, forcefield=forcefields["uff"], probe_symbol="He"):
    """
    Mix framework and probe UFF parameters with Lorentz-Berthelot rules for each framework atom.

    :param ase_atoms: the atoms object from ASE
    :param forcefield: a dictionary of force field parameters for the atoms in the framework
    :param probe_symbol: the symbol of the probe atom, defaults to He (optional)
    :return: The mixed sigma and epsilon values for the interactions.
    """

    # * let 's ge tteh lorentz berthelot rules mixed force field params for the interactions
    import numpy as np
    import operator

    ff = forcefield
    frame_ff = np.vstack(
        operator.itemgetter(*ase_atoms.get_chemical_symbols())(ff)
    )  # framework force field
    probe_ff = ff[probe_symbol]  # proble
    lb_sigma = (frame_ff[:, 1] + probe_ff[1]) / 2  # mixed sigma
    lb_epsilon = np.sqrt(frame_ff[:, 0] * probe_ff[0])  # mixed epsilon
    lb_ff = np.vstack((lb_epsilon, lb_sigma)).T  # mixed sigma an d epsilon
    return lb_ff


def show_percolated_paths_as_df(paths, path_windows, path_lengths):
    """
    Build a sorted DataFrame of pore-graph paths, window diameters along each path, and hop counts.
    
    :param paths: A list of paths, where each entry represent teh index of a pocket/node on the pore graph
    :param path_windows: Path windows are the sizes of the windows encountered along the path
    :param path_lengths: The number of jumps or edges along the path
    :return: The function `show_percolated_paths_as_df` is returning a pandas DataFrame sorted by the
    path lengths, where the columns show the paths, path windows, and path lengths."""

    import pandas as pd

    # * make a dataframe with the paths and the window sizes and the path lengths
    df = pd.DataFrame(
        {"paths": paths, "path_windows": path_windows, "path_lengths": path_lengths}
    )
    # Sort the DataFrame by the number of elements in the list in column 0, and then by column 2
    df_sorted = (
        df.assign(col0_len=df["paths"].apply(len))
        .sort_values(["col0_len", "path_lengths"])
        .drop("col0_len", axis=1)
    )

    # Print the sorted DataFrame
    # print(df_sorted)
    return df_sorted
