import requests
import time
import random
import re
from collections import defaultdict
import pickle
import equilibrator_api
import equilibrator_cache
import numpy as np
import json
import re
import equilibrator_api
from equilibrator_api import ComponentContribution
from tqdm import tqdm

# Step 1: Get all reaction IDs
reaction_list_url = "https://rest.kegg.jp/list/reaction"
response = requests.get(reaction_list_url)
reaction_ids = [line.split("\t")[0] for line in response.text.strip().split("\n")]

def convert_modelseed_equation(eq):
    """
    Convert ModelSEED/KEGG equation to equilibrator format.
    Example:
    "(1) cpd00001[0] + (1) cpd00012[0] <=> (2) cpd00009[0] + (1) cpd00067[0]"
    -> "1 cpd00001 + 1 cpd00012 <=> 2 cpd00009 + 1 cpd00067"
    """
    def clean_side(side):
        # Matches: (coeff) cpd[compartment]
        comps = []
        for part in side.split('+'):
            m = re.match(r"\s*\((\d+)\)\s*([a-zA-Z0-9_]+)\[\d+\]", part.strip())
            if m:
                coeff, cid = m.groups()
                comps.append(f"{coeff} {cid}")
        return ' + '.join(comps)
    if '<=>' in eq:
        left, right = eq.split('<=>')
    elif '=>' in eq:
        left, right = eq.split('=>')
    elif '<=' in eq:
        left, right = eq.split('<=')
    elif '>' in eq:
        left, right = eq.split('>')
    elif '<' in eq:
        left, right = eq.split('<')
    else:
        raise ValueError("Equation must contain '<=>' or '=>'")
    left_clean = clean_side(left)
    right_clean = clean_side(right)
    return f"{left_clean} <=> {right_clean}"

BALANCING_COMPOUNDS = [
    "C00001",  # H2O
    "C00080",  # H+
    "C00009",  # Pi
    "C00011",  # CO2
    "C00010",  # CoA
    "C00002",  # ATP   
    "C00008",  # ADP
    "C00020",  # AMP
    "C00003",  # NAD+
    "C00004",  # NADH
    "C00006",  # NADP+
    "C00005",  # NADPH
    "C00016",  # FAD
    "C01352",  # FADH2
    "C00035",  # GDP
    "C00044",  # GTP
]

CC = ComponentContribution()
def recalc_deltag(eq:str, CC:ComponentContribution, convert_from_modelseed:bool=False, verbose:bool=False):
    if convert_from_modelseed:
        eq = convert_modelseed_equation(eq)
    #- parse reaction
    rea = CC.parse_reaction_formula(eq)
    if not rea.is_balanced():
        if verbose: print(f"Warning: reaction not balanced: {eq}; attemting to balance")
        #- try to balance with BALANCING_COMPOUNDS:
        for cpd in BALANCING_COMPOUNDS:
            #print(cpd)
            try:
                rea.balance_with_compound(CC.get_compound(cpd))
            except Exception as e:
                continue # Silently skip to the next compound
            if rea.is_balanced():
                if verbose: print(f"Balanced by {cpd}")
                break
        if not rea.is_balanced():
            reab = CC.balance_by_oxidation(rea)
            if reab.is_balanced():
                if verbose: print(f"Balanced by oxidation")
                rea = reab
    if not rea.is_balanced():
        if verbose: print(f"Reaction still unbalanced.")
    res = CC.dg_prime(rea)
    return res.magnitude.nominal_value, res.magnitude.std_dev


#- parse ModelSEED reactions
with open('../downloads/reactions.json') as f:
    jdat = json.load(f)
    
MAP_KEGG_TO_DELTAG = {}
for entry in tqdm(jdat,desc="Processing ModelSEED reactions"):
    if entry['status'] in ['EMPTY', 'CPDFORMERROR']: #- skip deleted/invalid entries
        ##print(f"Skipping {entry['id']} with status {entry['status']}")
        continue
    eid = entry['id']
    kid = None
    ##print(eid)
    if entry['aliases'] is None:
        continue
    kid = None
    deltag = None
    for alias in entry['aliases']:
        m = re.match(r"KEGG:\s*(\S+)", alias)
        if m:
            kid = m.group(1)
            break
    deltag = entry['deltag']
    deltag_std = entry['deltagerr']
    if deltag is None or deltag == 10000000:
        ##print (f"Recalculating deltag for {eid} ({kid})")
        eq = entry['equation']
        if eq is not None:
            try:
                deltag, deltag_std = recalc_deltag(eq, CC, convert_from_modelseed=True)
                ##print(f"New deltag: {deltag} +/- {deltag_std}")
            except Exception as e:
                ##print(f"Error recalculating deltag for {eid} ({kid}): {e}")
                deltag, deltag_std = None, None
    if kid is not None and deltag is not None:
        if not kid in MAP_KEGG_TO_DELTAG:
            MAP_KEGG_TO_DELTAG[kid] = {'ids': [], 'deltags': [], 'detags_sdev': []}
        MAP_KEGG_TO_DELTAG[kid]['ids'].append(eid)
    if deltag is not None and kid is not None:
        MAP_KEGG_TO_DELTAG[kid]['deltags'].append(deltag)
        if deltag_std is not None:
            MAP_KEGG_TO_DELTAG[kid]['detags_sdev'].append(deltag_std)

#- save to file
with open("../res/kegg_to_modelseed_deltag.pkl", "wb") as f:
    pickle.dump(MAP_KEGG_TO_DELTAG, f)
    
##- rxn02733 error in recalculation
## rxn03385 formula confusion

#- iterate through MAP_KEGG_TO_MODELSEED and average deltags values (they usually are the same)
for kid in MAP_KEGG_TO_DELTAG:
    MAP_KEGG_TO_DELTAG[kid]['deltag'] = np.mean(MAP_KEGG_TO_DELTAG[kid]['deltags'])   
    

# Function to fetch reaction details in batches
def fetch_reaction_data(batch, max_retries=3, sleep_time=2):
    """Fetch reaction data from KEGG API, with error handling and retries."""
    url = f"https://rest.kegg.jp/get/{'+'.join(batch)}"
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=10)  # Set timeout to avoid hanging
            if response.status_code == 200 and response.text.strip():
                return response.text  # Return valid data
            print(
                f"Warning: Empty or failed response (attempt {attempt+1}/{max_retries})"
            )
        except requests.RequestException as e:
            print(f"Error: {e} (attempt {attempt+1}/{max_retries})")
        time.sleep(sleep_time)  # Wait before retrying
    print("Error: Failed to fetch reaction data after multiple attempts.")
    return None  # Return None to indicate failure


def parse_equation(equation):
    def extract_compounds(compound_str):
        """Extract compounds and their counts, keeping symbolic coefficients like 'n' or 'n+1'."""
        compound_list = []
        components = re.split(r" \+ ", compound_str)

        for comp in components:
            match = re.match(r"([\d\w\+\-\(\)\s]*)\s*([CG]\d+)", comp.strip())
            # Match numeric counts, symbolic expressions (e.g., n, n+1), or no count
            if match:
                count = match.group(1).strip()
                compound = match.group(2)

                if count == "":
                    count = "1"  # Assume 1 if no explicit coefficient

                compound_list.append((compound, count))  # Keep count as a string

        return compound_list

    # Split into reactants and products
    if " <=> " in equation:
        reactants, products = equation.split(" <=> ")
    elif " => " in equation:
        reactants, products = equation.split(" => ")
    else:
        return None

    return {"ipt": extract_compounds(reactants), "opt": extract_compounds(products)}


# Step 2: Process reactions in batches
batch_size = 16  # Adjust based on API limits
reaction_data = {}
failed_batches = []
for i in range(0, len(reaction_ids), batch_size):
    batch = reaction_ids[i : i + batch_size]

    data = fetch_reaction_data(batch)

    # - if data is None, append batch to failed_batches
    if data is None:
        failed_batches.append(batch)
        print(
            f"Failed to fetch batch {i//batch_size + 1} of {len(reaction_ids)//batch_size + 1}..."
        )
        continue

    # Parse equations
    entries = data.split("ENTRY       ")
    for entry in entries[1:]:  # First split is empty
        lines = entry.split("\n")
        rid = lines[0].split()[0]
        for line in lines:
            if line.startswith("EQUATION"):
                equation = line.replace("EQUATION    ", "").strip()
                parsed_eq = parse_equation(equation)
                if parsed_eq:
                    reaction_data[rid] = parsed_eq
                break
    # - every 10th iteration pickele reaction_data to disk
    if (i // batch_size + 1) % 10 == 0:
        with open("reaction_data.pkl", "wb") as f:
            pickle.dump(reaction_data, f)
            print(
                f"Fetched batch {i//batch_size + 1} of {len(reaction_ids)//batch_size + 1}..."
            )

    # - draw a random number between 0.1 and 0.75
    time.sleep(0.1 + 0.65 * random.random())  # Avoid rate limiting

# - write the final reaction_data to disk
with open("reaction_data.pkl", "wb") as f:
    pickle.dump(reaction_data, f)

# - write the failed_batches to disk
if len(failed_batches) > 0:
    with open("reactions_failed_batches.pkl", "wb") as f:
        pickle.dump(failed_batches, f)

# - report how many b atches failed
print(f"{len(failed_batches)} batches failed.")

print("done.")
