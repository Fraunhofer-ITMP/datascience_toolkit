# -*- coding: utf-8 -*-

import os
import pickle
import logging
import json
import pandas as pd
import networkx as nx
import requests
from collections import defaultdict

import streamlit as st
from stqdm import stqdm
from pybel import BELGraph
from pybel.dsl import Protein, Abundance, BiologicalProcess

from chembl_webresource_client.new_client import new_client

import plotly.graph_objects as go
from plotly.subplots import make_subplots


logger = logging.getLogger("__name__")

DATA_DIR = "data/"


def load_kg(path):
    infile = open(path, "rb")
    kg = pickle.load(infile)
    infile.close()
    return kg


def disease_figures(disease_name):
    """Function to generate figures for the disease overview."""

    disease_name = disease_name.lower()

    data_files = os.listdir(DATA_DIR + f"{disease_name}_kgg")
    graph_path = [file for file in data_files if file.endswith(".pkl")][0]

    # Basic graph stats
    graph = load_kg(DATA_DIR + f"{disease_name}_kgg/" + graph_path)

    nodes_data = {}
    for node in graph.nodes():
        if type(node).__name__ not in nodes_data:
            nodes_data[type(node).__name__] = 0

        nodes_data[type(node).__name__] += 1

    edge_dict = {}
    for source, target in graph.edges(data=False):
        s = type(source).__name__
        t = type(target).__name__

        if f"{s}-{t}" not in edge_dict and f"{t}-{s}" not in edge_dict:
            edge_dict[f"{s}-{t}"] = 0

        if f"{s}-{t}" in edge_dict:
            edge_dict[f"{s}-{t}"] += 1
        if f"{t}-{s}" in edge_dict:
            edge_dict[f"{t}-{s}"] += 1

    node_df = pd.DataFrame(nodes_data.items(), columns=["Node type", "Count"])
    node_df = node_df.sort_values(by="Count", ascending=False)

    edge_df = pd.DataFrame(edge_dict.items(), columns=["Edge type", "Count"])
    edge_df = edge_df.sort_values(by="Count", ascending=False)

    fig = make_subplots(rows=1, cols=2, subplot_titles=("Node summary", "Edge summary"))

    fig.append_trace(
        go.Bar(
            x=node_df["Node type"],
            y=node_df["Count"],
            text=node_df["Count"],
        ),
        row=1,
        col=1,
    )

    fig.append_trace(
        go.Bar(
            x=edge_df["Edge type"],
            y=edge_df["Count"],
            text=edge_df["Count"],
        ),
        row=1,
        col=2,
    )

    fig.update_layout(title="Graph summary", showlegend=False)

    # Update xaxis properties
    fig.update_xaxes(title_text="Node types", row=1, col=1)
    fig.update_xaxes(title_text="Edge types", row=1, col=2)

    # Update yaxis properties
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=2)

    st.plotly_chart(fig, use_container_width=True)

    # Namespace summary
    namespace_data = {}
    for node in graph.nodes():
        if node.namespace not in namespace_data:
            namespace_data[node.namespace] = 0

        namespace_data[node.namespace] += 1

    namespace_df = pd.DataFrame(namespace_data.items(), columns=["Namespace", "Count"])
    namespace_df = namespace_df.sort_values(by="Count", ascending=False)

    figure_ns = go.Figure()
    figure_ns.add_trace(
        go.Bar(
            x=namespace_df["Namespace"],
            y=namespace_df["Count"],
            text=namespace_df["Count"],
            marker_color="#ed872d",
        )
    )
    figure_ns.update_layout(title="Namespace summary")
    figure_ns.update_xaxes(title_text="Namespace")
    figure_ns.update_yaxes(title_text="Count")

    st.plotly_chart(figure_ns, use_container_width=True)

    # TODO: Add drug specific information


def RetMech(chemblIds) -> list:
    """Function to retrieve mechanism of actions and target proteins from ChEMBL

    :param chemblIds:
    :return:
    """
    # Create a placeholder for displaying current progress
    getMech = new_client.mechanism

    mechList = []
    for chemblid in chemblIds:
        mechs = getMech.filter(molecule_chembl_id=chemblid).only(
            ["mechanism_of_action", "target_chembl_id", "action_type"]
        )

        mechList.append(list(mechs))

    named_mechList = dict(zip(chemblIds, mechList))
    named_mechList = {k: v for k, v in named_mechList.items() if v}
    return named_mechList


def RetAct(chemblIds) -> dict:
    """Function to retrieve associated assays from ChEMBL

    :param chemblIds:
    :return:
    """
    GetAct = new_client.activity
    getTar = new_client.target
    ActList = []
    filtered_list = [
        "assay_chembl_id",
        "assay_type",
        "pchembl_value",
        "target_chembl_id",
        "target_organism",
        "bao_label",
        "target_type",
    ]

    for chembl in chemblIds:
        acts = GetAct.filter(
            molecule_chembl_id=chembl,
            pchembl_value__isnull=False,
            assay_type_iregex="(B|F)",
            target_organism="Homo sapiens",
        ).only(filtered_list)

        data = []

        for d in acts:

            if float(d.get("pchembl_value")) < 6:
                continue

            if d.get("bao_label") != "single protein format":
                continue

            tar = d.get("target_chembl_id")
            tar_dict = getTar.get(tar)

            try:
                if tar_dict["target_type"] in ("CELL-LINE", "UNCHECKED"):
                    continue
            except KeyError:
                continue

            data.append(d)

        ActList.append(list(data))

    named_ActList = dict(zip(chemblIds, ActList))
    named_ActList = {k: v for k, v in named_ActList.items() if v}
    return named_ActList


def RetDrugInd(chemblIDs) -> dict:
    """Function to retrieve associated diseases from ChEMBL

    :param chemblIDs:
    :return:
    """
    getDrugInd = new_client.drug_indication

    drugIndList = []
    for chemblid in chemblIDs:
        drugInd = getDrugInd.filter(molecule_chembl_id=chemblid).only("mesh_heading")
        drugIndList.append(list(drugInd))

    named_drugIndList = dict(zip(chemblIDs, drugIndList))
    named_drugIndList = {k: v for k, v in named_drugIndList.items() if v}
    return named_drugIndList


def Ret_chembl_protein(sourceList) -> list:
    """Method to retrieve ChEMBL ids which are proteins/targets

    :param sourceList:
    :return:
    """
    protein_List = []
    for item in sourceList:
        for j in range(len(sourceList[item])):
            protein_List.append(sourceList[item][j]["target_chembl_id"])

    protein_List = set(protein_List)
    protein_List = list(filter(None, protein_List))
    return protein_List


def chembl2uniprot(chemblIDs) -> dict:
    """Method to convert ChEMBL id to UNIPROT and get associated REACTOME pathways

    :param chemblIDs:
    :return:
    """
    getTarget = new_client.target
    chem2Gene2path = []
    chemHasNoPath = set()
    chemNotprotein = set()

    chem2path = defaultdict(list)

    # Loop to ensure it is a protein
    for chemblid in chemblIDs:
        chem = getTarget.filter(chembl_id=chemblid).only("target_components")

        try:
            uprot_id = chem[0]["target_components"][0]["accession"]

            if not uprot_id:
                chemHasNoPath.add(chemblid)

        except IndexError:
            chemHasNoPath.add(chemblid)

    logger.info(f"No UniProt information available for {len(chemHasNoPath)} proteins.")

    chemblIDs_filtered = [item for item in chemblIDs if item not in chemHasNoPath]

    # Get gene symbol from ChEMBL and filtering the list for human proteins only
    for chemblid in chemblIDs_filtered:

        chem = getTarget.filter(chembl_id=chemblid).only("target_components")
        getGene = chem[0]["target_components"][0]["target_component_synonyms"]
        try:
            getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]

            if not getGene:
                chemNotprotein.add(chemblid)

        except IndexError:
            chemNotprotein.add(chemblid)

    chemblIDs_filtered = [
        item for item in chemblIDs_filtered if item not in chemNotprotein
    ]

    # Extracting data for valid proteins only
    for chemblid in chemblIDs_filtered:
        chem = getTarget.filter(chembl_id=chemblid).only("target_components")

        # UniProt data
        uprot_id = chem[0]["target_components"][0]["accession"]

        # Gene symbol
        getGene = chem[0]["target_components"][0]["target_component_synonyms"]
        getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]

        # Pathway data
        chem2path = [
            item
            for item in chem[0]["target_components"][0]["target_component_xrefs"]
            if item["xref_src_db"] == "Reactome"
        ]

        uprot = {"accession": uprot_id}
        chem2path.append(uprot)
        chem2path.append(getGene)
        chem2Gene2path.append(chem2path)

    named_chem2Gene2path = dict(zip(chemblIDs_filtered, chem2Gene2path))
    named_chem2Gene2path = {k: v for k, v in named_chem2Gene2path.items() if v}
    return named_chem2Gene2path


def chembl2gene2path(chem2geneList, ActList):
    """Method for updating chembl protein nodes with gene symbol.

    :param chem2geneList:
    :param ActList:
    :return:
    """
    for item in chem2geneList:
        sizeOfitem = len(chem2geneList[item])
        gene = chem2geneList[item][sizeOfitem - 1]["component_synonym"]
        for jtem in ActList:
            for i in range(len(ActList[jtem])):
                if item == ActList.get(jtem)[i]["target_chembl_id"]:
                    newkey = {"Protein": gene}
                    ActList[jtem][i].update(newkey)

    return ActList


def chem2moa_rel(named_mechList, org, graph: BELGraph) -> BELGraph:
    """Method to create the graph"""
    pos = ["POSITIVE ALLOSTERIC MODULATOR", "AGONIST", "ACTIVATOR", "PARTIAL AGONIST"]
    neg = ["INHIBITOR", "NEGATIVE ALLOSTERIC MODULATOR", "ANTAGONIST", "BLOCKER"]
    misc = [
        "MODULATOR",
        "DISRUPTING AGENT",
        "SUBSTRATE",
        "OPENER",
        "SEQUESTERING AGENT",
    ]

    for chembl_name, chembl_entries in named_mechList.items():
        for info in chembl_entries:
            graph.add_association(
                Abundance(namespace="ChEMBL", name=chembl_name),
                BiologicalProcess(namespace="MOA", name=info["mechanism_of_action"]),
                citation="ChEMBL database",
                evidence="ChEMBL query",
            )

            if not info["target_chembl_id"]:
                continue

            if "Protein" in info:
                if info["action_type"] in pos:
                    graph.add_increases(
                        Abundance(namespace="ChEMBL", name=chembl_name),
                        Protein(namespace=org, name=info["Protein"]),
                        citation="ChEMBL database",
                        evidence="ChEMBL query",
                    )
                if info["action_type"] in neg:
                    graph.add_decreases(
                        Abundance(namespace="ChEMBL", name=chembl_name),
                        Protein(namespace=org, name=info["Protein"]),
                        citation="ChEMBL database",
                        evidence="ChEMBL query",
                    )

                if info["action_type"] in misc:
                    graph.add_association(
                        Abundance(namespace="ChEMBL", name=chembl_name),
                        Protein(namespace=org, name=info["Protein"]),
                        citation="ChEMBL database",
                        evidence="ChEMBL query",
                    )

    return graph


def chem2act_rel(named_ActList, org, graph: BELGraph) -> BELGraph:
    """Method to add bioassay edges to the KG."""
    for chemical, chem_entries in named_ActList.items():
        for chem_data in chem_entries:
            if chem_data["target_chembl_id"]:
                if "Protein" in chem_data:
                    graph.add_association(
                        Abundance(
                            namespace="ChEMBLAssay", name=chem_data["assay_chembl_id"]
                        ),
                        Protein(namespace=org, name=chem_data["Protein"]),
                        citation="ChEMBL database",
                        evidence="ChEMBL query",
                    )

            graph.add_association(
                Abundance(namespace="ChEMBL", name=chemical),
                Abundance(namespace="ChEMBLAssay", name=chem_data["assay_chembl_id"]),
                citation="ChEMBL database",
                evidence="ChEMBL query",
                annotation={
                    "assayType": chem_data["assay_type"],
                    "pChEMBL": chem_data["pchembl_value"],
                },
            )

    return graph


def gene2path_rel(named_chem2geneList, org, graph) -> BELGraph:
    """Method to add protein and reactome data to KG"""
    for item in named_chem2geneList:
        itemLen = len(named_chem2geneList[item]) - 1
        for j in range(itemLen - 1):
            graph.add_association(
                Protein(
                    namespace=org,
                    name=named_chem2geneList[item][itemLen]["component_synonym"],
                ),
                BiologicalProcess(
                    namespace="Reactome", name=named_chem2geneList[item][j]["xref_name"]
                ),
                citation="ChEMBL database",
                evidence="ChEMBL query",
                annotation={
                    "Reactome": "https://reactome.org/content/detail/"
                    + named_chem2geneList[item][j]["xref_id"]
                },
            )

    return graph


def uniprot_rel(named_uprotList, org, graph) -> BELGraph:
    """Method to add UniProt related edges"""
    for item in named_uprotList:
        fun = list(named_uprotList[item]["Function"].keys())
        bp = list(named_uprotList[item]["BioProcess"].keys())
        for f in fun:
            if str(named_uprotList[item]["Gene"]) != "nan" and not isinstance(
                named_uprotList[item]["Gene"], dict
            ):
                graph.add_association(
                    Protein(namespace=org, name=named_uprotList[item]["Gene"]),
                    BiologicalProcess(namespace="GOMF", name=f),
                    citation="UniProt database",
                    evidence="UniProt query",
                )
            else:
                graph.add_association(
                    Protein(namespace=org, name=item),
                    BiologicalProcess(namespace="GOMF", name=f),
                    citation="UniProt database",
                    evidence="UniProt query",
                )

        for b in bp:
            if str(named_uprotList[item]["Gene"]) != "nan" and not isinstance(
                named_uprotList[item]["Gene"], dict
            ):
                graph.add_association(
                    Protein(namespace=org, name=named_uprotList[item]["Gene"]),
                    BiologicalProcess(namespace="GOBP", name=b),
                    citation="UniProt database",
                    evidence="UniProt query",
                )
            else:
                graph.add_association(
                    Protein(namespace=org, name=item),
                    BiologicalProcess(namespace="GOBP", name=b),
                    citation="UniProt database",
                    evidence="UniProt query",
                )

        if str(named_uprotList[item]["Gene"]) != "nan" and not isinstance(
            named_uprotList[item]["Gene"], dict
        ):
            nx.set_node_attributes(
                graph,
                {
                    Protein(
                        namespace=org, name=named_uprotList[item]["Gene"]
                    ): "https://3dbionotes.cnb.csic.es/?queryId="
                    + item
                },
                "3Dbio",
            )

            nx.set_node_attributes(
                graph,
                {
                    Protein(
                        namespace=org, name=named_uprotList[item]["Gene"]
                    ): "https://www.uniprot.org/uniprotkb/"
                    + item
                },
                "UniProt",
            )

        else:
            nx.set_node_attributes(
                graph,
                {
                    Protein(
                        namespace=org, name=item
                    ): "https://3dbionotes.cnb.csic.es/?queryId="
                    + item
                },
                "3Dbio",
            )

            nx.set_node_attributes(
                graph,
                {
                    Protein(
                        namespace=org, name=item
                    ): "https://www.uniprot.org/uniprotkb/"
                    + item
                },
                "UniProt",
            )

    return graph


def searchDisease(keyword):
    """Finding disease identifiers using OpenTargets API"""
    disease_name = str(keyword)

    query_string = """
        query searchAnything ($disname:String!){
            search(queryString:$disname,entityNames:"disease",page:{size:20,index:0}){
                total
                hits {
                    id
                    entity
                    name
                    description
                }
            }
        }
        """

    variables = {"disname": disease_name}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    df = pd.DataFrame(api_response["data"]["search"]["hits"])
    df = df[df["entity"] == "disease"]
    df.drop(columns=["entity"], inplace=True)

    return df


def getDrugCount(disease_id):
    """Finding the number of drugs associated with a disease using OpenTargets API"""
    efo_id = disease_id

    query_string = """
    query associatedTargets($my_efo_id: String!){
        disease(efoId: $my_efo_id){
            id
            name
            knownDrugs{
                uniqueTargets
                uniqueDrugs
                count
            }
        }
    }

    """

    # Set variables object of arguments to be passed to endpoint
    variables = {"my_efo_id": efo_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    # get the count value from api_repsonse dict
    api_response = api_response["data"]["disease"]["knownDrugs"]["count"]
    return api_response


def GetDiseaseAssociatedDrugs(disease_id, CT_phase):
    """Finding drugs associated with a disease using OpenTargets API"""
    efo_id = disease_id
    size = getDrugCount(efo_id)

    query_string = """
    query associatedTargets($my_efo_id: String!, $my_size: Int){
        disease(efoId: $my_efo_id){
            id
            name
            knownDrugs(size:$my_size){
                uniqueTargets
                uniqueDrugs
                count
                rows{
                    approvedSymbol
                    approvedName
                    prefName
                    drugType
                    drugId
                    phase
                    ctIds
                }

            }
        }
    }
    """

    # Set variables object of arguments to be passed to endpoint
    variables = {"my_efo_id": efo_id, "my_size": size}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    df = pd.DataFrame(api_response["data"]["disease"]["knownDrugs"]["rows"])

    if df.empty:
        st.write("No drugs found in clinical trials")
        return df

    df = df.loc[df["phase"] >= int(CT_phase), :]
    df["id"] = efo_id
    df["disease"] = api_response["data"]["disease"]["name"]
    return df


def GetDiseaseAssociatedProteins(disease_id):
    """Finding proteins associated with a disease using OpenTargets API"""
    efo_id = str(disease_id)

    query_string = """
        query associatedTargets{
            disease(efoId: $efo_id){
                id
                name
                associatedTargets(page:{size:15000,index:0}){
                    count
                    rows {
                        target {
                            id
                            approvedSymbol
                            proteinIds {
                                id
                                source
                            }
                        }
                        score
                    }
                }
            }
        }
    """

    # replace $efo_id with value from efo_id
    query_string = query_string.replace("$efo_id", f'"{efo_id}"')

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string})

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    temp_list = []
    for item in api_response["data"]["disease"]["associatedTargets"]["rows"]:
        for obj in item["target"]["proteinIds"]:
            if obj["source"] == "uniprot_swissprot":
                # st.write(obj)
                uprot = obj["id"]
                source = obj["source"]
                score = item["score"]
                ensg = item["target"]["id"]
                name = item["target"]["approvedSymbol"]
                temp = {
                    "Protein": name,
                    "ENSG": ensg,
                    "UniProt": uprot,
                    "Source": source,
                    "Score": score,
                }
                temp_list.append(temp)

    df = pd.DataFrame(temp_list)
    df["disease_id"] = efo_id

    return df


def GetDiseaseAssociatedProteinsPlot(df):
    """Plotting the protein confidence scores associated with a disease."""
    st.markdown("**Protein-Disease Association summary**")
    st.markdown(
        f"""We have identified {len(df)} proteins (Swiss-Prot) associated with the disease. Please note that the proteins identified may not be unique if you combined two or more diseases. Following is a histogram that shows distribution of proteins based on scores provided by OpenTargets. The scores are influenced by various factors such as genetic associations, expression, mutations, known pathways, targeting drugs and so on."""
    )

    prot_fig = go.Figure()
    prot_fig.add_trace(
        go.Bar(
            x=df["Protein"].head(20),
            y=df["Score"].head(20),
            marker=dict(color="#a4d3b3"),
            name="Protein",
        )
    )
    prot_fig.update_layout(
        title="Distribution of top 20 proteins based on OpenTargets score",
        xaxis_title="Protein",
        yaxis_title="Score",
    )
    st.plotly_chart(prot_fig, use_container_width=True)

    score = st.number_input(
        "Enter threshold score (recommended > 0.3):",
        min_value=0.0,
        max_value=1.0,
        value=0.3,
        step=0.1,
    )
    return score


def ExtractFromUniProt(uniprot_id) -> dict:
    """Uniprot parser to retrieve information about OMIM disease, reactome pathway, biological process,
     and molecular functions.

    :param uniprot_id:
    :return:
    """
    Uniprot_Dict = []

    mapped_uprot = []

    for id in stqdm(uniprot_id, "Extracting data from UniProt"):
        # Retrieve data for id in text format if found in uniprot
        ret_uprot = requests.get(
            "https://www.uniprot.org/uniprot/" + id + ".txt"
        ).text.split("\n")

        if ret_uprot == [""]:
            continue

        id_copy = id
        mapped_uprot.append(id_copy)

        k = 0
        id = {}
        id["Disease"] = {}
        id["Reactome"] = {}
        id["Function"] = {}
        id["BioProcess"] = {}
        id["Gene"] = {}

        # parse each line looking for info about disease, pathway, funcn, bp and so on
        for line in ret_uprot:
            # parse lines with disease and extract disease names and omim ids
            if "-!- DISEASE:" in line:
                if "[MIM:" in line:
                    dis = line.split(":")
                    id["Disease"].update({dis[1][1:-5]: dis[2][:-1]})

            # extract reactome ids and names
            if "Reactome;" in line:
                ract = line.split(";")
                id["Reactome"].update({ract[2][1:-2]: ract[1][1:]})

            # look for functions
            if " F:" in line:
                fn = line.split(";")
                id["Function"].update({fn[2][3:]: fn[1][1:]})

            # look for biological processes
            if " P:" in line and "GO;" in line:
                bp = line.split(";")
                id["BioProcess"].update({bp[2][3:]: bp[1][1:]})

            if "GN   Name" in line:
                if k == 0:
                    gene = line.split("=")
                    gene = gene[1].split(" ")
                    if ";" in gene[0]:
                        gene = gene[0].split(";")
                        gene = {"Gene": gene[0]}
                    else:
                        gene = {"Gene": gene[0]}
                    id.update(gene)
                    k += 1

        Uniprot_Dict.append(id)

    Uniprot_Dict = dict(zip(mapped_uprot, Uniprot_Dict))

    return Uniprot_Dict


def createInitialKG(session_inputs):
    """Creating the initial Knowledge Graph using the disease and protein data."""
    efo_id = session_inputs["disease_id"]
    ct_phase = session_inputs["ct_phase"]

    for functions in stqdm(
        ["disease_drugs", "disease_proteins", "disease_snp"],
        "Fetching real-time data from databases. Be patient!",
    ):
        if functions == "disease_drugs":
            drugs_df = GetDiseaseAssociatedDrugs(efo_id, ct_phase)
            drugs_df = drugs_df.reset_index(drop=True)

        elif functions == "disease_proteins":
            dis2prot_df = GetDiseaseAssociatedProteins(efo_id)
            dis2prot_df = dis2prot_df.reset_index(drop=True)

    return drugs_df, dis2prot_df


def finalizeKG(filtered_protein_df: pd.DataFrame, session_inputs: dict):
    """Finalizing the Knowledge Graph by adding proteins and drugs data."""

    # create empty KG
    kg = BELGraph(name=session_inputs["kg_name"], version="0.0.1")

    for metadata_functions in "prot", "cmpds":
        if metadata_functions == "prot":
            unique_proteins = list(set(filtered_protein_df["UniProt"]))
            uprot_ext = ExtractFromUniProt(unique_proteins)
        elif metadata_functions == "cmpds":
            drugs_df = st.session_state.drugs_df
            if not drugs_df.empty:

                for rel_function_1 in stqdm(
                    ["chembl2mech", "chembl2act"], desc="Fetching chemical relations"
                ):
                    if rel_function_1 == "chembl2mech":
                        chembl2mech = RetMech(list(set(drugs_df["drugId"])))
                    elif rel_function_1 == "chembl2act":
                        chembl2act = RetAct(list(set(drugs_df["drugId"])))

                prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(
                    chembl2mech
                )
                prtn_as_chembl = set(prtn_as_chembl)
                chembl2uprot = chembl2uniprot(prtn_as_chembl)

                for rel_function_2 in stqdm(
                    ["chembl2act", "chembl2mech"], desc="Fetching protein relations"
                ):
                    if rel_function_2 == "chembl2act":
                        chembl2act = chembl2gene2path(chembl2uprot, chembl2act)
                    elif rel_function_2 == "chembl2mech":
                        chembl2mech = chembl2gene2path(chembl2uprot, chembl2mech)

                kg = chem2moa_rel(chembl2mech, "HGNC", kg)
                kg = chem2act_rel(chembl2act, "HGNC", kg)
                kg = gene2path_rel(chembl2uprot, "HGNC", kg)
                kg = uniprot_rel(uprot_ext, "HGNC", kg)

    st.write("Your KG is now generated!", "\n")
    return kg
