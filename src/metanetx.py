"""MetaNetX RDF parsing and graph construction."""

import gc
import hashlib
import pickle
import re
from collections import defaultdict
from pathlib import Path

import networkx as nx
from rdflib import Graph, Namespace
from rdflib.namespace import RDFS
from tqdm import tqdm

MNX = Namespace("https://rdf.metanetx.org/schema/")
KEGG_RXN_PREFIX = "https://identifiers.org/kegg.reaction:"
KEGG_C_PREFIXES = (
    "https://identifiers.org/kegg.compound:",
    "https://www.genome.jp/kegg/compound/",
)
RE_KEGG_C = re.compile(r"keggC:(C\d+)")

CACHE_VERSION = "v1"
CACHE_DIR_NAME = ".cache"
CACHE_FILE_TEMPLATE = "{stem}.{ttl_hash}.indices.pkl"
MD5_CHUNK_SIZE = 64 * 1024 * 1024


def _compute_md5(file_path: Path, chunk_size: int = MD5_CHUNK_SIZE) -> str:
    """Compute an MD5 hash for a file using chunked reads."""
    hasher = hashlib.md5()
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _cache_path(ttl_path: Path, ttl_hash: str) -> Path:
    """Return the cache path for a MetaNetX TTL hash."""
    cache_dir = ttl_path.parent / CACHE_DIR_NAME
    cache_dir.mkdir(parents=True, exist_ok=True)
    filename = CACHE_FILE_TEMPLATE.format(stem=ttl_path.stem, ttl_hash=ttl_hash)
    return cache_dir / filename


def _load_cached_indices(cache_path: Path, ttl_hash: str):
    """Load cached indices if they match the requested TTL hash."""
    if not cache_path.exists():
        return None

    try:
        with cache_path.open("rb") as handle:
            cached = pickle.load(handle)
    except Exception as exc:  # pragma: no cover - cache corruption path
        print(f"Warning: Failed to load MetaNetX cache {cache_path}: {exc}")
        return None

    if cached.get("ttl_hash") != ttl_hash or cached.get("cache_version") != CACHE_VERSION:
        return None

    return cached


def _store_indices_cache(cache_path: Path, indices: dict) -> None:
    """Persist MetaNetX indices for reuse."""
    try:
        with cache_path.open("wb") as handle:
            pickle.dump(indices, handle, protocol=pickle.HIGHEST_PROTOCOL)
    except OSError as exc:  # pragma: no cover - filesystem error path
        print(f"Warning: Failed to write MetaNetX cache {cache_path}: {exc}")


def load_metanetx_indices(ttl_path: Path) -> dict:
    """
    Parse MetaNetX TTL file and build lookup indices (with caching).

    Returns dict with keys:
        - kegg_to_mnxr: {kegg_rxn_id: {mnxr_ids}}
        - mnx_chem_to_kegg: {mnx_chem_uri: {kegg_compound_ids}}
        - rxn_comments: {kegg_rxn_id: comment_text}
        - part_data: {part_uri: {chem, coef, comp}}
        - rxn_structure: {mnxr_id: {left: [parts], right: [parts]}}
        - ttl_hash: MD5 hash of the TTL file used to generate the cache
        - cache_version: cache metadata version
    """
    ttl_hash = _compute_md5(ttl_path)
    cache_path = _cache_path(ttl_path, ttl_hash)

    cached = _load_cached_indices(cache_path, ttl_hash)
    if cached is not None:
        print(f"Using cached MetaNetX indices: {cache_path}")
        return cached

    print(f"Loading Turtle file: {ttl_path}")
    g = Graph()
    g.parse(str(ttl_path), format="turtle")

    print("Indexing graph data...")
    kegg_to_mnxr = defaultdict(set)
    mnx_chem_to_kegg = defaultdict(set)
    rxn_comments = {}
    part_data = defaultdict(dict)
    rxn_structure = defaultdict(lambda: {"left": [], "right": []})

    for s, p, o in tqdm(g, desc="Scanning Triples"):
        s_str, o_str = str(s), str(o)

        if p == MNX.reacXref:
            if o_str.startswith(KEGG_RXN_PREFIX):
                kegg_id = o_str.rsplit(":", 1)[-1]
                mnx_id = s_str.rsplit("/", 1)[-1]
                kegg_to_mnxr[kegg_id].add(mnx_id)

        elif p == MNX.chemXref:
            for prefix in KEGG_C_PREFIXES:
                if o_str.startswith(prefix):
                    kegg_cid = o_str.rsplit(":", 1)[-1]
                    mnx_chem_to_kegg[s].add(kegg_cid)
                    break

        elif p == RDFS.comment:
            if s_str.startswith(KEGG_RXN_PREFIX):
                kegg_id = s_str.rsplit(":", 1)[-1]
                rxn_comments[kegg_id] = str(o)

        elif p == MNX.left or p == MNX.right:
            mnx_id = s_str.rsplit("/", 1)[-1]
            side = "left" if p == MNX.left else "right"
            rxn_structure[mnx_id][side].append(o)

        elif p == MNX.chem:
            part_data[s]["chem"] = o
        elif p == MNX.coef:
            part_data[s]["coef"] = float(o)
        elif p == MNX.comp:
            part_data[s]["comp"] = str(o).rsplit("/", 1)[-1]

    print("Freeing RDF graph memory...")
    del g
    gc.collect()

    indices = {
        "kegg_to_mnxr": dict(kegg_to_mnxr),
        "mnx_chem_to_kegg": dict(mnx_chem_to_kegg),
        "rxn_comments": rxn_comments,
        "part_data": dict(part_data),
        "rxn_structure": {k: dict(v) for k, v in rxn_structure.items()},
    }

    indices["ttl_hash"] = ttl_hash
    indices["cache_version"] = CACHE_VERSION

    _store_indices_cache(cache_path, indices)

    return indices


def _chem_to_kegg_ids(chem_uri, comment_text: str, mnx_chem_to_kegg: dict) -> set:
    """Resolve MNX chemical URI to KEGG compound IDs."""
    candidates = mnx_chem_to_kegg.get(chem_uri, set())
    if len(candidates) <= 1 or not comment_text:
        return candidates
    tokens = set(RE_KEGG_C.findall(comment_text))
    preferred = candidates & tokens
    return preferred or candidates


def build_kegg_digraph(
    kegg_rxn_ids: set[str],
    indices: dict,
) -> nx.DiGraph:
    """
    Build a directed graph from KEGG reaction IDs using MetaNetX structure.

    Nodes have 'group' attribute: 'Reaction' or 'Compound'
    Edges have 'coef' and 'compartment' attributes.
    """
    kegg_to_mnxr = indices["kegg_to_mnxr"]
    rxn_structure = indices["rxn_structure"]
    part_data = indices["part_data"]
    rxn_comments = indices["rxn_comments"]
    mnx_chem_to_kegg = indices["mnx_chem_to_kegg"]

    G = nx.DiGraph()

    for kegg_rxn in tqdm(kegg_rxn_ids, desc="Building graph"):
        mnxrs = kegg_to_mnxr.get(kegg_rxn)
        if not mnxrs:
            continue

        rxn_node = f"kegg:{kegg_rxn}"
        if rxn_node not in G:
            G.add_node(rxn_node, group="Reaction", name=kegg_rxn)

        comment = rxn_comments.get(kegg_rxn, "")

        for mnxr in mnxrs:
            structure = rxn_structure.get(mnxr)
            if not structure:
                continue

            for side in ["left", "right"]:
                for part_uri in structure.get(side, []):
                    p_attrs = part_data.get(part_uri)
                    if not p_attrs or "chem" not in p_attrs:
                        continue

                    kegg_cids = _chem_to_kegg_ids(
                        p_attrs["chem"], comment, mnx_chem_to_kegg
                    )

                    for cid in kegg_cids:
                        cmp_node = f"kegg:{cid}"
                        if cmp_node not in G:
                            G.add_node(cmp_node, group="Compound", name=cid)

                        edge_attrs = {
                            "coef": p_attrs.get("coef", 1),
                            "compartment": p_attrs.get("comp"),
                        }

                        if side == "left":
                            G.add_edge(cmp_node, rxn_node, **edge_attrs)
                        else:
                            G.add_edge(rxn_node, cmp_node, **edge_attrs)

    return G
