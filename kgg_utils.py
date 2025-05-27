# -*- coding: utf-8 -*-

import io
import json
import logging
import os
import pickle
import tempfile
import zipfile
from collections import Counter, defaultdict

import networkx as nx
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
import pybel
import pybel.struct.mutation.induction as induction
import pybel_jupyter
import requests
import streamlit as st
from chembl_webresource_client.new_client import new_client
from pandasgwas import get_variants
from pandasgwas.get_variants import get_variants_by_efo_id
from plotly.subplots import make_subplots
from pybel import BELGraph
from pybel.dsl import Abundance, BiologicalProcess, Gene, Pathology, Protein
from rdkit import Chem
from rdkit.Chem import Descriptors
from stqdm import stqdm

logger = logging.getLogger("__name__")

DATA_DIR = "data/"


def load_kg(path):
    infile = open(path, "rb")
    kg = pickle.load(infile)
    infile.close()
    return kg


state = st.session_state


def disease_figures(disease_name, graph: BELGraph = None):
    """Function to generate figures for the disease overview."""

    if graph is None:
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

    if "figures" not in state:
        state["figures"] = {}

    state["figures"]["graph_summary"] = fig
    state["figures"]["namespace_summary"] = figure_ns

    # TODO: Add drug specific information


def RetMech(chemblIds) -> list:
    """Function to retrieve mechanism of actions and target proteins from ChEMBL

    :param chemblIds:
    :return:
    """
    # Create a placeholder for displaying current progress
    getMech = new_client.mechanism

    mechList = []
    for chemblid in stqdm(chemblIds, desc="Fetching mechanims of actions"):
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

    for chembl in stqdm(chemblIds, desc="Fetching targets and associated assays"):
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
        # uniprot id is located in 2nd last pos
        uprot = chem2geneList[item][-2]["accession"]
        for jtem in ActList:
            for i in range(len(ActList[jtem])):
                if item == ActList.get(jtem)[i]["target_chembl_id"]:
                    newkey = {"Protein": gene, "Accession": uprot}
                    ActList[jtem][i].update(newkey)

    return ActList


# def chembl2gene2path(chem2geneList, ActList):
#     """Method for updating chembl protein nodes with gene symbol.
#
#     :param chem2geneList:
#     :param ActList:
#     :return:
#     """
#     for item in chem2geneList:
#         sizeOfitem = len(chem2geneList[item])
#         gene = chem2geneList[item][sizeOfitem - 1]["component_synonym"]
#         for jtem in ActList:
#             for i in range(len(ActList[jtem])):
#                 if item == ActList.get(jtem)[i]["target_chembl_id"]:
#                     newkey = {"Protein": gene}
#                     ActList[jtem][i].update(newkey)
#
#     return ActList


def chem2moa_rel(named_mechList, org, otp_prots, graph: BELGraph) -> BELGraph:
    """Method to create the monkeypox graph

    :param named_mechList:
    :param org:
    :param graph: BEL graph of Monkeypox
    :return:
    """

    # identified types of chemical and protein action types
    # ['INHIBITOR','NEGATIVE ALLOSTERIC MODULATOR','POSITIVE ALLOSTERIC MODULATOR','ANTAGONIST','AGONIST','MODULATOR','BLOCKER','ACTIVATOR','DISRUPTING AGENT', 'SUBSTRATE', 'OPENER','PARTIAL AGONIST','SEQUESTERING AGENT']
    # following lists are used to determine type of edge relationships
    pos = ["POSITIVE ALLOSTERIC MODULATOR", "AGONIST", "ACTIVATOR", "PARTIAL AGONIST"]
    neg = ["INHIBITOR", "NEGATIVE ALLOSTERIC MODULATOR", "ANTAGONIST", "BLOCKER"]
    misc = [
        "MODULATOR",
        "DISRUPTING AGENT",
        "SUBSTRATE",
        "OPENER",
        "SEQUESTERING AGENT",
    ]

    for chembl_name, chembl_entries in stqdm(
        named_mechList.items(), desc="Populating Chemical-MoA edges"
    ):
        for info in chembl_entries:
            graph.add_qualified_edge(
                Abundance(namespace="ChEMBL", name=chembl_name),
                BiologicalProcess(namespace="MOA", name=info["mechanism_of_action"]),
                relation="hasMechanismOfAction",
                citation="ChEMBL database",
                evidence="ChEMBL query",
            )

            if not info["target_chembl_id"]:
                continue

            if (
                "Protein" in info
                and "Accession" in info
                and info["Accession"] in otp_prots
            ):
                if info["action_type"] in pos and info["Accession"] in otp_prots:
                    graph.add_increases(
                        Abundance(namespace="ChEMBL", name=chembl_name),
                        Protein(namespace=org, name=info["Protein"]),
                        citation="ChEMBL database",
                        evidence="ChEMBL query",
                    )
                if info["action_type"] in neg and info["Accession"] in otp_prots:
                    graph.add_decreases(
                        Abundance(namespace="ChEMBL", name=chembl_name),
                        Protein(namespace=org, name=info["Protein"]),
                        citation="ChEMBL database",
                        evidence="ChEMBL query",
                    )

                if info["action_type"] in misc and info["Accession"] in otp_prots:
                    graph.add_qualified_edge(
                        Abundance(namespace="ChEMBL", name=chembl_name),
                        Protein(namespace=org, name=info["Protein"]),
                        relation="targets",
                        citation="ChEMBL database",
                        evidence="ChEMBL query",
                    )

            # else:
            # graph.add_association(
            # Abundance(namespace='ChEMBL', name=chembl_name),
            # Protein(namespace=org, name=info['target_chembl_id']),
            # citation='ChEMBL database',
            # evidence='ChEMBL query'
            # )

    return graph


def chem2act_rel(named_ActList, org, otp_prots, graph: BELGraph) -> BELGraph:
    """Method to add bioassay edges to the KG.

    :param named_ActList:
    :param org:
    :param graph:
    :return:
    """
    for chemical, chem_entries in stqdm(
        named_ActList.items(), desc="Adding bioassay edges to BEL"
    ):
        for chem_data in chem_entries:
            if chem_data["target_chembl_id"]:
                if (
                    "Protein" in chem_data
                    and "Accession" in chem_data
                    and chem_data["Accession"] in otp_prots
                ):
                    graph.add_qualified_edge(
                        Abundance(
                            namespace="ChEMBLAssay", name=chem_data["assay_chembl_id"]
                        ),
                        Protein(namespace=org, name=chem_data["Protein"]),
                        relation="hasTarget",
                        citation="ChEMBL database",
                        evidence="ChEMBL query",
                    )
                # else:
                # graph.add_association(
                # Abundance(namespace='ChEMBLAssay', name=chem_data['assay_chembl_id']),
                # Protein(namespace=org, name=chem_data['target_chembl_id']),
                # citation='ChEMBL database',
                # evidence='ChEMBL query'
                # )

            graph.add_qualified_edge(
                Abundance(namespace="ChEMBL", name=chemical),
                Abundance(namespace="ChEMBLAssay", name=chem_data["assay_chembl_id"]),
                relation="hasAssay",
                citation="ChEMBL database",
                evidence="ChEMBL query",
                annotation={
                    "assayType": chem_data["assay_type"],
                    "pChEMBL": chem_data["pchembl_value"],
                },
            )

    return graph


def gene2path_rel(named_chem2geneList, org, otp_prots, graph) -> BELGraph:
    """Method to add protein and reactome data to KG

    :param named_chem2geneList:
    :param org:
    :param graph:
    :return:
    """
    for item in named_chem2geneList:
        itemLen = len(named_chem2geneList[item]) - 1
        for j in range(itemLen - 1):
            # checks if uprot id is in otp_proteins
            if named_chem2geneList[item][itemLen - 1]["accession"] in otp_prots:
                graph.add_qualified_edge(
                    Protein(
                        namespace=org,
                        name=named_chem2geneList[item][itemLen]["component_synonym"],
                    ),
                    BiologicalProcess(
                        namespace="Reactome",
                        name=named_chem2geneList[item][j]["xref_name"],
                    ),
                    relation="hasPathway",
                    citation="ChEMBL database",
                    evidence="ChEMBL query",
                    annotation={
                        "Reactome": "https://reactome.org/content/detail/"
                        + named_chem2geneList[item][j]["xref_id"]
                    },
                )

    return graph


def uniprot_rel(named_uprotList, org, graph) -> BELGraph:
    """Method to add UniProt related edges

    :param named_uprotList:
    :param org:
    :param graph:
    :return:
    """
    for item in stqdm(named_uprotList, desc="Populating Uniprot edges"):
        fun = list(named_uprotList[item]["Function"].keys())
        bp = list(named_uprotList[item]["BioProcess"].keys())
        for f in fun:
            if str(named_uprotList[item]["Gene"]) != "nan" and not isinstance(
                named_uprotList[item]["Gene"], dict
            ):
                graph.add_qualified_edge(
                    Protein(namespace=org, name=named_uprotList[item]["Gene"]),
                    BiologicalProcess(namespace="GOMF", name=f),
                    relation="hasMolecularFunction",
                    citation="UniProt database",
                    evidence="UniProt query",
                )
            else:
                graph.add_qualified_edge(
                    Protein(namespace=org, name=item),
                    BiologicalProcess(namespace="GOMF", name=f),
                    relation="hasMolecularFunction",
                    citation="UniProt database",
                    evidence="UniProt query",
                )

        for b in bp:
            if str(named_uprotList[item]["Gene"]) != "nan" and not isinstance(
                named_uprotList[item]["Gene"], dict
            ):
                graph.add_qualified_edge(
                    Protein(namespace=org, name=named_uprotList[item]["Gene"]),
                    BiologicalProcess(namespace="GOBP", name=b),
                    relation="hasBiologicalProcess",
                    citation="UniProt database",
                    evidence="UniProt query",
                )
            else:
                graph.add_qualified_edge(
                    Protein(namespace=org, name=item),
                    BiologicalProcess(namespace="GOBP", name=b),
                    relation="hasBiologicalProcess",
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
                    ): "https://3dbionotes.cnb.csic.es/?queryId=" + item
                },
                "3Dbio",
            )

            nx.set_node_attributes(
                graph,
                {
                    Protein(
                        namespace=org, name=named_uprotList[item]["Gene"]
                    ): "https://www.uniprot.org/uniprotkb/" + item
                },
                "UniProt",
            )

        else:
            nx.set_node_attributes(
                graph,
                {
                    Protein(
                        namespace=org, name=item
                    ): "https://3dbionotes.cnb.csic.es/?queryId=" + item
                },
                "3Dbio",
            )

            nx.set_node_attributes(
                graph,
                {
                    Protein(
                        namespace=org, name=item
                    ): "https://www.uniprot.org/uniprotkb/" + item
                },
                "UniProt",
            )

    return graph


def searchDisease(keyword):
    """Finding disease identifiers using OpenTargets API"""
    disease_name = str(keyword)
    if disease_name is None:
        st.error("Please resubmit the disease.")
        st.stop()
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
    # if r is None:
    #     st.error("Please resubmit the disease.")
    #     st.stop()
    api_response = json.loads(r.text)

    df = pd.DataFrame(api_response["data"]["search"]["hits"])
    if not df.empty:
        df = df[df["entity"] == "disease"]
        df.drop(columns=["entity"], inplace=True)
        return df
    else:
        st.warning("Disease with given keyword is not found. Please try again!")
        st.stop()


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
    efo_id = state.get("disease_id", "")
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


@st.cache_data(ttl=3600, show_spinner="Fetching protein data...")
def GetDiseaseAssociatedProteins(efo_id):
    """Fixed version that properly handles cache invalidation"""
    #    st.write(f"🔄 Disease Lookup: {efo_id}")

    if not efo_id:
        st.error("No disease ID provided")
        return pd.DataFrame()

    base_url = "https://api.platform.opentargets.org/api/v4/graphql"
    query = """
    query associatedTargets($efoId: String!, $index: Int!) {
      disease(efoId: $efoId) {
        associatedTargets(page: {size: 3000, index: $index}) {
          rows {
            target {
              id
              approvedSymbol
              proteinIds { id source }
            }
            score
          }
        }
      }
    }
    """

    all_results = []
    index = 0

    while True:
        response = requests.post(
            base_url,
            json={"query": query, "variables": {"efoId": efo_id, "index": index}},
        )
        response.raise_for_status()
        data = response.json()

        batch = data["data"]["disease"]["associatedTargets"]["rows"]
        if not batch:
            break

        all_results.extend(batch)
        index += 1

    processed = []
    for item in all_results:
        target = item["target"]
        for protein in target["proteinIds"]:
            if protein["source"] == "uniprot_swissprot":
                processed.append(
                    {
                        "Protein": target["approvedSymbol"],
                        "ENSG": target["id"],
                        "UniProt": protein["id"],
                        "Source": protein["source"],
                        "Score": item["score"],
                        "disease_id": efo_id,
                    }
                )

    df = pd.DataFrame(processed)
    st.write(f"✅ Retrieved {len(df)} proteins for {efo_id}")
    return df


# def GetDiseaseAssociatedProteins(disease_id):
# """Finding proteins associated with a disease using OpenTargets API"""
# efo_id = str(disease_id)

# query_string = """
# query associatedTargets{
# disease(efoId: $efo_id){
# id
# name
# associatedTargets(page:{size:15000,index:0}){
# count
# rows {
# target {
# id
# approvedSymbol
# proteinIds {
# id
# source
# }
# }
# score
# }
# }
# }
# }
# """

# # replace $efo_id with value from efo_id
# query_string = query_string.replace("$efo_id", f'"{efo_id}"')

# # Set base URL of GraphQL API endpoint
# base_url = "https://api.platform.opentargets.org/api/v4/graphql"

# # Perform POST request and check status code of response
# r = requests.post(base_url, json={"query": query_string})

# # Transform API response from JSON into Python dictionary and print in console
# api_response = json.loads(r.text)

# temp_list = []
# for item in api_response["data"]["disease"]["associatedTargets"]["rows"]:
# for obj in item["target"]["proteinIds"]:
# if obj["source"] == "uniprot_swissprot":
# # st.write(obj)
# uprot = obj["id"]
# source = obj["source"]
# score = item["score"]
# ensg = item["target"]["id"]
# name = item["target"]["approvedSymbol"]
# temp = {
# "Protein": name,
# "ENSG": ensg,
# "UniProt": uprot,
# "Source": source,
# "Score": score,
# }
# temp_list.append(temp)

# df = pd.DataFrame(temp_list)
# df["disease_id"] = efo_id

# return df


def GetDiseaseAssociatedProteinsPlot(df):
    """Plotting the protein confidence scores associated with a disease."""
    if "figures" in st.session_state and "protein_score" in st.session_state.figures:
        del st.session_state.figures["protein_score"]

    st.markdown("**Protein-Disease Association summary**")
    st.markdown(
        f"""We have identified {len(df)} proteins (Swiss-Prot) associated with the disease. 
        Please note that the proteins identified may not be unique if you combined two or more diseases. 
        Following is a histogram that shows distribution of proteins based on scores provided by OpenTargets. 
        The scores are influenced by various factors such as genetic associations, expression, 
        mutations, known pathways, targeting drugs and so on."""
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

    current_disease = st.session_state.get("user_disease", "Unknown Disease")
    #    st.write(f"The disease name is {current_disease}")

    prot_fig.update_layout(
        title=f"Distribution of top 20 proteins for {current_disease} based on OpenTargets score",
        xaxis_title="Protein",
        yaxis_title="Score",
    )

    st.plotly_chart(prot_fig, use_container_width=True)
    st.caption(f"Top 20 proteins based on OpenTargets score")

    if "figures" not in st.session_state:
        st.session_state.figures = {}
    st.session_state.figures["protein_score"] = prot_fig


def clearPlotsandGraphs():
    """
    This function clears the plots and graphs from the session state.
    """
    if "figures" in state:
        del state["figures"]

    if "graphs" in state:
        del state["graphs"]

    if "graph_summary" in state:
        del state["graph_summary"]

    if "namespace_summary" in state:
        del state["namespace_summary"]


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


def GetDiseaseSNPs():
    try:
        disease_id = state.get("disease_id", "")
        snps = get_variants_by_efo_id(disease_id)
        snps_df = snps.genomic_contexts
        snps_df["disease_id"] = disease_id
        # snps_functional_class_df = snps.variants
        # snps_functional_class_df['disease_id'] = disease_id
        snps_df = snps_df.reset_index(drop=True)
        return snps_df

    except:
        st.write(
            "No SNPs found. This could be either because 1. The identifier of your disease of interest is not compatible with Experimental Factor Ontology (EFO) or 2. No SNPs have been reported for the disease."
        )
        # snps_df = pd.DataFrame()
        # return(snps_df)


def GetViralProteins(query_disease):
    # file downloaded from https://www.genome.jp/ftp/db/virushostdb Dated: 12/09/2023
    virus = pd.read_csv(
        "https://raw.githubusercontent.com/Fraunhofer-ITMP/kgg/main/data/misc/virushostdb.csv"
    )
    # virus = pd.read_csv('../data/misc/virushostdb.csv')

    cols = ["virus tax id", "virus name", "DISEASE", "host tax id"]
    virus = virus[cols]

    # print(virus)

    # filter virus with host humans
    virus = virus.loc[virus["host tax id"] == 9606.0, :]
    virus = virus.reset_index(drop=True)

    # replace 9606.0 to 9606
    virus["host tax id"] = pd.to_numeric(virus["host tax id"], downcast="integer")

    # get the initial keyword for disease search
    # disease = GetQuery()

    # subset df with disease keyword
    virus_subset_1 = virus[
        virus["DISEASE"].str.contains(query_disease, na=False, case=False)
    ]

    if not virus_subset_1.empty:
        # print(disease)
        # print('\n')
        st.write(
            "The workflow has identified your query as a viral disease. In the next steps we will add viral proteins (SWISS-Prot) in the KG.",
            "\n",
        )

        # time.sleep(0.1)
        virus_name = st.text_input(
            "Would you like to look further for a specific virus? Please type its name or skip by typing 'no'.",
            placeholder="no",
        )

        # virus_name = st.selectbox(
        #     "Would you like to look further for a specific virus?",
        #     ("Yes", "No"),
        #     index=None,
        #     placeholder="Select Yes or No",
        # )

        st.write("You selected:", virus_name)

        # virus_name = input(
        #     'Do you want to look further for a specific virus? Please type its name or skip it by typing \'no\': ')

        # subset df with virus name
        if virus_name.lower() != "no":
            virus_subset_2 = virus[
                virus["virus name"].str.contains(virus_name, na=False, case=False)
            ]
            # print(virus_subset_2)
            # break
        else:
            virus_subset_2 = pd.DataFrame()

        # merge subsets of df_1 and df_2
        virus_subset_merge = pd.concat([virus_subset_1, virus_subset_2])

        virus_subset_merge = virus_subset_merge.drop_duplicates(keep="first")

        virus_subset_merge = virus_subset_merge.reset_index(drop=True)

        virus_subset_merge["index"] = virus_subset_merge.index

        virus_subset_merge.style.hide(axis="index")

        virus_subset_merge = virus_subset_merge[
            ["index", "virus tax id", "virus name", "DISEASE", "host tax id"]
        ]
        virus_subset_merge["virus tax id"] = virus_subset_merge["virus tax id"].astype(
            str
        )

        st.dataframe(virus_subset_merge)
        # display(HTML(virus_subset_merge.to_html(index=False)))

        temp_id = st.text_input(
            "Enter the index value(s). If multiple, use space, for example -> 0 1 3: ",
            placeholder=0,
        )

        # time.sleep(0.1)
        # temp_id = input('Enter the index value(s). If multiple, use space, for example -> 0 1 3: ')

        # print('\n')

        temp_id = temp_id.split(" ")
        temp_id = [int(x) for x in temp_id if x.strip()]
        # print(virus_subset_merge.loc[0]['virus tax id'])

        uprot_list = []

        for item in temp_id:
            tax_id = virus_subset_merge.loc[item]["virus tax id"]
            # print(tax_id)

            # fetch tax id related proteins from Uniprot
            # the link can be created from downloads option in uniprot
            query_string = (
                "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Clength%2Cgene_primary%2Cprotein_name&format=tsv&query=%28%28taxonomy_id%3A"
                + str(tax_id)
                + "%29+AND+%28reviewed%3Atrue%29%29"
            )

            # query_string = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Clength%2Cgene_primary%2Cprotein_name&format=tsv&query=%28%28taxonomy_id%3A11676%29%29+AND+%28reviewed%3Atrue%29'

            query_uniprot = requests.get(query_string)
            query_uniprot = query_uniprot.text.split("\n")

            query_uniprot_df = pd.DataFrame(
                [x.strip().split("\t") for x in query_uniprot]
            )
            cols = query_uniprot_df.iloc[0]
            # print(cols)
            query_uniprot_df = query_uniprot_df[1 : len(query_uniprot_df) - 1]
            query_uniprot_df.columns = cols
            temp = list(query_uniprot_df["Entry"])
            # print(len(temp))
            uprot_list.append(temp)

        uprot_list = [item for sublist in uprot_list for item in sublist]

        # st.write('A total of', str(len(uprot_list)), 'viral proteins have been identified.', '\n')
        st.write(
            f"A total of {str(len(uprot_list))} viral proteins have been identified."
        )

        return uprot_list


def getProtfromKG(mainGraph):
    prot_list = []
    for u, v, data in stqdm(
        mainGraph.edges(data=True), desc="Filtering Proteins/Genes"
    ):
        if "HGNC" in u.namespace:
            if u.name not in prot_list:
                prot_list.append(u.name)

        if "HGNC" in v.namespace:
            if v.name not in prot_list:
                prot_list.append(v.name)

    return prot_list


def snp2gene_rel(snp_df, graph):
    # print(snp_df.head(2))

    kg_prots = getProtfromKG(graph)

    unique_prots_df = pd.DataFrame(kg_prots, columns=["Proteins"])

    snp_df = unique_prots_df.merge(
        snp_df, how="inner", left_on="Proteins", right_on="gene.geneName"
    )

    # take SNPs which are within the gene sequence
    snp_df = snp_df.loc[snp_df["distance"] == 0]
    snp_df = snp_df.reset_index(drop=True)

    print(
        "A total of "
        + str(len(snp_df))
        + " SNPs have been identified from GWAS Central. Now adding relevant data"
    )
    print("\n")

    # graph = pybel.BELGraph(name='test', version="0.0.1")

    for i in stqdm(range(len(snp_df)), desc="Adding disease associated SNPs"):
        graph.add_qualified_edge(
            Protein(namespace="HGNC", name=snp_df["Proteins"][i]),
            Gene(namespace="dbSNP", name=snp_df["rsId"][i]),
            relation="hasGeneticVariant",
            citation="GWAS Central",
            evidence="SNPs for queried disease",
        )

    return graph


@st.cache_data(ttl=5, show_spinner="Fetching new disease data...")
def createInitialKG(_ct_phase):
    """Creating the initial Knowledge Graph using the disease and protein data."""
    efo_id = state.get("disease_id", "")
    if not efo_id:
        st.error("No disease ID found in session state.")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    ct_phase = _ct_phase
    drugs_df = pd.DataFrame()
    dis2prot_df = pd.DataFrame()
    dis2snp = pd.DataFrame()
    #    st.write(f"Disease ID inside createinitialkg: {efo_id}")
    for functions in stqdm(
        ["disease_drugs", "disease_proteins", "disease_snp"],
        "Fetching real-time data from databases. Be patient!",
    ):
        if functions == "disease_drugs":
            st.write("Fetching Drugs")
            drugs_df = GetDiseaseAssociatedDrugs(efo_id, ct_phase)
            drugs_df = drugs_df.reset_index(drop=True)

        elif functions == "disease_proteins":
            st.write("Fetching Proteins")
            dis2prot_df = GetDiseaseAssociatedProteins(state.get("disease_id", ""))
            dis2prot_df = dis2prot_df.reset_index(drop=True)

        elif functions == "disease_snp":
            st.write("Fetching SNPs")
            dis2snp = GetDiseaseSNPs()
            # dis2snp = dis2snp.reset_index(drop=True)

    return drugs_df, dis2prot_df, dis2snp


def getAdverseEffectCount(chembl_id):
    get_id = chembl_id

    query_string = """
        query AdverseEventsQuery(
          $chemblId: String!
          $index: Int = 0
          $size: Int = 10
        ) {
          drug(chemblId: $chemblId) {
            id
            maxLlr: adverseEvents(page: { index: 0, size: 1 }) {
              rows {
                logLR
              }
            }
            adverseEvents(page: { index: $index, size: $size }) {
              criticalValue
              count
              rows {
                name
                count
                logLR
                meddraCode
              }
            }
          }
        }

    """

    # Set variables object of arguments to be passed to endpoint
    variables = {"chemblId": get_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    # get the count value from api_repsonse dict
    api_response = api_response["data"]["drug"]["adverseEvents"]["count"]
    return api_response


def GetAdverseEvents(chem_list):
    api_response = pd.DataFrame()

    for chem in stqdm(chem_list, desc="Retrieving Adverse Effects for each drug"):
        chembl_id = chem

        try:
            # get total no. of adverse effects for a given drug
            count = getAdverseEffectCount(chembl_id)

            query_string = """
                query AdverseEventsQuery(
                  $chemblId: String!
                  $index: Int = 0
                  $size: Int!
                ) {
                  drug(chemblId: $chemblId) {
                    id
                    maxLlr: adverseEvents(page: { index: 0, size: 1 }) {
                      rows {
                        logLR
                      }
                    }
                    adverseEvents(page: { index: $index, size: $size }) {
                      criticalValue
                      count
                      rows {
                        name
                        count
                        logLR
                        meddraCode
                      }
                    }
                  }
                }

        """
            # Set variables object of arguments to be passed to endpoint
            variables = {"chemblId": chembl_id, "size": count}

            # Set base URL of GraphQL API endpoint
            base_url = "https://api.platform.opentargets.org/api/v4/graphql"

            # Perform POST request and check status code of response
            r = requests.post(
                base_url, json={"query": query_string, "variables": variables}
            )
            # r = requests.post(base_url, json={"query": query_string})
            # print(r.status_code)

            # Transform API response from JSON into Python dictionary and print in console
            api_response_temp = json.loads(r.text)

            api_response_temp = api_response_temp["data"]["drug"]["adverseEvents"][
                "rows"
            ]
            api_response_temp = pd.DataFrame(api_response_temp)
            api_response_temp["chembl_id"] = chembl_id

            api_response = pd.concat([api_response, api_response_temp])

        except:
            continue

    api_response.reset_index(drop=True, inplace=True)
    return api_response


def chembl2adverseEffect_rel(chembl_adveff_df, graph: BELGraph) -> BELGraph:
    """

    :param chembl_adveff_df:
    :param graph:
    :return:
    """

    for i in range(len(chembl_adveff_df)):
        graph.add_qualified_edge(
            Abundance(namespace="ChEMBL", name=str(chembl_adveff_df["chembl_id"][i])),
            Pathology(namespace="SideEffect", name=str(chembl_adveff_df["name"][i])),
            relation="hasAdverseEffect",
            citation="OpenTargets Platform",
            evidence="DrugReactions",
        )

    return graph


def getNodeList(nodeName, graph):
    # import pybel
    node_list = []
    for node in graph.nodes():
        if isinstance(node, pybel.dsl.Abundance):
            if node.namespace == nodeName:
                node_list.append(node.name)
    return node_list


def chembl_annotation(graph):
    chemblids = getNodeList("ChEMBL", graph)
    for item in stqdm(chemblids, desc="Adding ChEMBL URLs"):
        nx.set_node_attributes(
            graph,
            {
                Abundance(
                    namespace="ChEMBL", name=item
                ): "https://www.ebi.ac.uk/chembl/compound_report_card/" + item
            },
            "ChEMBL",
        )
    return graph


def chembl_name_annotation(graph, drugs_df):
    # create a dict {'Metformin':'CHEMBL1431'} for adding node attributes
    drugName_chembl_dict = {}
    for v, k in drugs_df[["prefName", "drugId"]].values:
        # print(k,v)
        drugName_chembl_dict.update({k: v})

    molecule = new_client.molecule

    # get all chemicals/drugs from kg
    chemblids = getNodeList("ChEMBL", graph)

    for chem in stqdm(chemblids, desc="Adding preferred names and trade names"):
        # print(chem)

        trade_names = []

        # fetch prefName and synonym from chembl
        getNames = molecule.filter(molecule_chembl_id=chem).only(
            ["pref_name", "molecule_synonyms"]
        )

        # get preferred name #not preferred because it can be empty sometimes
        # pref_name = getNames[0]['pref_name']

        # get preferred name of a drug from openTarget diseaseAssociatedDrug file #always has a preferredName
        nx.set_node_attributes(
            graph,
            {Abundance(namespace="ChEMBL", name=chem): drugName_chembl_dict[chem]},
            "PreferredName",
        )

        # get trade names and append to a list
        try:
            for item in getNames[0]["molecule_synonyms"]:
                if item["syn_type"] == "TRADE_NAME":
                    trade_names.append(item["molecule_synonym"])

            nx.set_node_attributes(
                graph,
                {Abundance(namespace="ChEMBL", name=chem): trade_names},
                "TradeName",
            )

        except:
            continue

    return graph


def getGeneOntolgyNodes(nodeName, graph):
    # import pybel
    node_list = []
    for node in graph.nodes():
        if isinstance(node, pybel.dsl.BiologicalProcess):
            # print(node)
            if node.namespace == nodeName:
                node_list.append(node.name)
    return node_list


def gene_ontology_annotation(graph, uprotDict):
    gobp_dict = {}
    gomf_dict = {}

    # create a merged dict of {bioprocess: bp_id} for mapping
    for prot in uprotDict:
        gobp_dict.update(uprotDict[prot]["BioProcess"])

    for prot in uprotDict:
        gomf_dict.update(uprotDict[prot]["Function"])

    # extract bps and mfs from kg
    bp_kg = getGeneOntolgyNodes("GOBP", graph)
    mf_kg = getGeneOntolgyNodes("GOMF", graph)

    for item in stqdm(bp_kg, desc="adding biological process annotations"):
        gobp_id = gobp_dict[item]
        # print(gobp_id)
        # print('https://www.ebi.ac.uk/QuickGO/term/'+gobp_id)
        nx.set_node_attributes(
            graph,
            {
                BiologicalProcess(
                    namespace="GOBP", name=item
                ): "https://www.ebi.ac.uk/QuickGO/term/" + gobp_id
            },
            "QuickGO",
        )
        nx.set_node_attributes(
            graph,
            {BiologicalProcess(namespace="GOBP", name=item): gobp_id},
            "Gene Ontology identifier",
        )

    for item in stqdm(mf_kg, desc="adding molecular function annotations"):
        gomf_id = gomf_dict[item]
        # print(gobp_id)
        # print('https://www.ebi.ac.uk/QuickGO/term/'+gobp_id)
        nx.set_node_attributes(
            graph,
            {
                BiologicalProcess(
                    namespace="GOMF", name=item
                ): "https://www.ebi.ac.uk/QuickGO/term/" + gomf_id
            },
            "QuickGO",
        )
        nx.set_node_attributes(
            graph,
            {BiologicalProcess(namespace="GOMF", name=item): gomf_id},
            "Gene Ontology identifier",
        )

    return graph


def finalizeKG(filtered_protein_df: pd.DataFrame, session_inputs: dict):
    """Finalizing the Knowledge Graph by adding proteins and drugs data."""

    # create empty KG
    kg = BELGraph(name=session_inputs["kg_name"], version="0.0.1")

    for metadata_functions in "prot", "cmpds", "snps":
        if metadata_functions == "prot":
            unique_proteins = list(set(filtered_protein_df["UniProt"]))
            # st.write(unique_proteins)
            uprot_ext = ExtractFromUniProt(unique_proteins)
            kg = uniprot_rel(uprot_ext, "HGNC", kg)
            kg = gene_ontology_annotation(kg, uprot_ext)

            if "human_protein" not in state:
                state["human_protein"] = uprot_ext
            elif uprot_ext != state["human_protein"]:
                state["human_protein"] = uprot_ext

            viral_prot = state.viral_prot
            if viral_prot:
                vir_uprot_ext = ExtractFromUniProt(viral_prot)
                # st.write(vir_uprot_ext)
                kg = uniprot_rel(vir_uprot_ext, "VP", kg)

                # kg = gene_ontology_annotation(kg, vir_uprot_ext)

                if "viral_protein" not in state:
                    state["viral_protein"] = vir_uprot_ext
                elif vir_uprot_ext != state["viral_protein"]:
                    state["viral_protein"] = vir_uprot_ext

        elif metadata_functions == "cmpds":
            drugs_df = state.drugs_df
            if not drugs_df.empty:
                # for rel_function_1 in stqdm(
                #     ["chembl2mech", "chembl2act"], desc="Fetching chemical relations"
                # ):
                #     if rel_function_1 == "chembl2mech":
                #         chembl2mech = RetMech(list(set(drugs_df["drugId"])))
                #     elif rel_function_1 == "chembl2act":
                #         chembl2act = RetAct(list(set(drugs_df["drugId"])))
                #
                # prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(chembl2mech)
                # prtn_as_chembl = set(prtn_as_chembl)
                # chembl2uprot = chembl2uniprot(prtn_as_chembl)
                #
                # for rel_function_2 in stqdm(
                #     ["chembl2act", "chembl2mech"], desc="Fetching protein relations"
                # ):
                #     if rel_function_2 == "chembl2act":
                #         chembl2act = chembl2gene2path(chembl2uprot, chembl2act)
                #     elif rel_function_2 == "chembl2mech":
                #         chembl2mech = chembl2gene2path(chembl2uprot, chembl2mech)
                #         st.write('Here')
                #         st.write(chembl2mech)

                st.write(
                    "A total of "
                    + str(len(list(set(drugs_df["drugId"]))))
                    + " drugs have been identified. Now fetching relevant data"
                )

                chembl2mech = RetMech(list(set(drugs_df["drugId"])))
                # st.write("Mech",chembl2mech)
                chembl2act = RetAct(list(set(drugs_df["drugId"])))
                # st.write('Act',chembl2act)

                prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(
                    chembl2mech
                )
                prtn_as_chembl = set(prtn_as_chembl)
                prtn_as_chembl = list(prtn_as_chembl)
                chembl2uprot = chembl2uniprot(prtn_as_chembl)
                # st.write('Chem2prot',chembl2uprot)

                chembl2act_2 = chembl2gene2path(chembl2uprot, chembl2act)
                chembl2mech_2 = chembl2gene2path(chembl2uprot, chembl2mech)

                # st.write("Mech2", chembl2act_2)
                # st.write('Act2', chembl2act_2)

                kg = chem2moa_rel(chembl2mech_2, "HGNC", unique_proteins, kg)
                kg = chem2act_rel(chembl2act_2, "HGNC", unique_proteins, kg)
                kg = gene2path_rel(chembl2uprot, "HGNC", unique_proteins, kg)

                adv_effect = GetAdverseEvents(list(set(drugs_df["drugId"])))

                # viral_prot = kgg_utils.GetViralProteins(state["user_disease"])
                #
                # if "viral_prot" not in state:
                #     state["viral_prot"] = viral_prot
                # elif not state['viral_prot'] == viral_prot:
                #     state["viral_prot"] = viral_prot

                if "adv_effect" not in state:
                    state["adv_effect"] = adv_effect
                elif not state["adv_effect"].equals(adv_effect):
                    state["adv_effect"] = adv_effect

                kg = chembl2adverseEffect_rel(adv_effect, kg)

                kg = chembl_name_annotation(kg, drugs_df)

        elif metadata_functions == "snps":
            dis2snp_df = state.dis2snp_df
            # st.write(state)
            # if not dis2snp_df.empty:
            #     kg = snp2gene_rel(dis2snp_df, kg)
            try:
                kg = snp2gene_rel(dis2snp_df, kg)
            except:
                continue

    #    st.write(state)

    st.write("Your KG is now generated!", "\n")
    return kg


def get_graph_summary(graph):
    """Printing summary similar to PyBEL summary."""
    rv_basic = [
        ("Name", graph.name),
        ("Version", graph.version),
    ]

    rv_stats = [
        ("Nodes", graph.number_of_nodes()),
        ("Namespaces", len(graph.count.namespaces())),
        ("Edges", graph.number_of_edges()),
        ("Annotations", len(graph.count.annotations())),
        ("Citations", graph.number_of_citations()),
        ("Authors", graph.number_of_authors()),
        ("Components", nx.number_weakly_connected_components(graph)),
        ("Warnings", graph.number_of_warnings()),
        ("Network Density", "{:.2E}".format(nx.density(graph))),
    ]

    return rv_basic, rv_stats


def GetSmiles(drugs_df, colname):
    # chembl ids can be antibodies which have sequence instead of smiles

    drugs = drugs_df
    # drugs = pd.read_csv('data/kgs/metabolic diseases/t1dm/diseaseAssociatedDrugs.csv')
    drugs_list = set(list(drugs[colname]))

    molecule = new_client.molecule
    temp_list = []
    unused_list = []

    for item in stqdm(
        drugs_list, "Getting SMILES for CHEMBL ids and generating descriptors"
    ):
        try:
            mol = molecule.filter(chembl_id=str(item)).only(
                ["molecule_chembl_id", "molecule_structures"]
            )
            mol_smiles = mol[0]["molecule_structures"]["canonical_smiles"]
            temp = [item, mol_smiles]
            temp_list.append(temp)

        except:
            # st.write(f'Id {str(item)} could not be parsed')
            # print(item)
            unused_list.append(item)
            continue

    temp_df = pd.DataFrame(temp_list, columns=["drugId", "smiles"])

    unused_df = pd.DataFrame(unused_list, columns=["Unparsed_drugs"])

    return (temp_df, unused_df)


def ro5_filter(df):
    temp_df = []

    violation_counts_list = []

    pass_list = []

    for item in df["smiles"]:
        molecule = Chem.MolFromSmiles(item)

        #         Lipinski Ro5

        #         Moleculer Weight <= 500
        #         LogP <= 5
        #         H-Bond Donor Count <= 5
        #         H-Bond Acceptor Count <= 10

        molecular_weight = Descriptors.MolWt(molecule)

        h_bond_acceptors = Descriptors.NOCount(molecule)

        h_bond_donor = Descriptors.NHOHCount(molecule)

        logp = Descriptors.MolLogP(molecule)

        conditions = [
            molecular_weight <= 500,
            h_bond_acceptors <= 10,
            h_bond_donor <= 5,
            logp <= 5,
        ]

        violation_counts = conditions.count(False)

        violation_counts_list.append(violation_counts)

        descriptors = [molecular_weight, h_bond_acceptors, h_bond_donor, logp]

        temp_df.append(descriptors)

        pass_ro5 = conditions.count(True) >= 3

        if pass_ro5 == True:
            pass_list.append(0)

        else:
            pass_list.append(1)

    names = "MW HBA HBD LogP".split(" ")
    temp_df = pd.DataFrame(temp_df, columns=names)

    df["Violation(s)_ro5"] = violation_counts_list
    df["Lipinski_ro5"] = pass_list

    df = pd.concat([df, temp_df], axis=1)

    return df


def ghose_filter(df):
    temp_df = []

    pass_list = []

    for item in df["smiles"]:
        molecule = Chem.MolFromSmiles(item)

        # ghose descriptors
        #     Molecular weight between 160 and 480
        #     LogP between -0.4 and +5.6
        #     Atom count between 20 and 70
        #     Molar refractivity between 40 and 130

        molar_refractivity = Chem.Crippen.MolMR(molecule)

        number_of_atoms = Chem.rdchem.Mol.GetNumAtoms(molecule)

        molecular_weight = Descriptors.MolWt(molecule)

        logp = Descriptors.MolLogP(molecule)

        conditions = [
            molecular_weight >= 160 and molecular_weight <= 480,
            number_of_atoms >= 20 and number_of_atoms <= 70,
            molar_refractivity >= 40 and molar_refractivity <= 130,
            logp >= -0.4 and logp <= 5.6,
        ]

        pass_ghose = conditions.count(True) == 4

        descriptors = [number_of_atoms, molar_refractivity]

        temp_df.append(descriptors)

        if pass_ghose == True:
            pass_list.append(0)

        else:
            pass_list.append(1)

    names = "AtomNum MolRefractivity".split(" ")
    temp_df = pd.DataFrame(temp_df, columns=names)

    df["Ghose"] = pass_list

    df = pd.concat([df, temp_df], axis=1)

    return df


def veber_filter(df):
    temp_df = []

    pass_list = []

    for item in df["smiles"]:
        molecule = Chem.MolFromSmiles(item)

        #         Veber filter
        #         Rotatable bonds <= 10
        #         Topological polar surface area <= 140

        topological_surface_area_mapping = Chem.QED.properties(molecule).PSA

        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)

        conditions = [rotatable_bonds <= 10, topological_surface_area_mapping <= 140]

        pass_veber = conditions.count(True) == 2

        descriptors = [rotatable_bonds, topological_surface_area_mapping]

        temp_df.append(descriptors)

        if pass_veber == True:
            pass_list.append(0)

        else:
            pass_list.append(1)

    names = "RotBond TPSA".split(" ")
    temp_df = pd.DataFrame(temp_df, columns=names)

    df["Veber"] = pass_list

    df = pd.concat([df, temp_df], axis=1)

    return df


def reos_filter(df):
    temp_df = []

    pass_list = []

    for item in df["smiles"]:
        molecule = Chem.MolFromSmiles(item)

        # REOS:
        #     Molecular weight between 200 and 500
        #     LogP between -5.0 and +5.0
        #     H-bond donor count between 0 and 5
        #     H-bond acceptor count between 0 and 10
        #     Formal charge between -2 and +2
        #     Rotatable bond count between 0 and 8
        #     Heavy atom count between 15 and 50

        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)
        h_bond_donor = Descriptors.NumHDonors(molecule)
        h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
        formal_charge = Chem.rdmolops.GetFormalCharge(molecule)
        heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(molecule)

        # print(molecular_weight,logp,h_bond_donor,h_bond_acceptors,formal_charge,rotatable_bonds,heavy_atoms)

        conditions = [
            molecular_weight >= 200 and molecular_weight <= 500,
            logp >= -0.4 and logp <= 5.6,
            h_bond_donor >= 0 and h_bond_donor <= 5,
            h_bond_acceptors >= 0 and h_bond_acceptors <= 10,
            formal_charge >= -2 and formal_charge <= 2,
            rotatable_bonds >= 0 and rotatable_bonds <= 8,
            heavy_atoms >= 15 and heavy_atoms <= 50,
        ]

        descriptors = [formal_charge, heavy_atoms]

        temp_df.append(descriptors)

        # print(conditions)
        pass_reos = conditions.count(True) == 7

        if pass_reos == True:
            pass_list.append(0)

        else:
            pass_list.append(1)

    names = "Charge HeavyAtom".split(" ")
    temp_df = pd.DataFrame(temp_df, columns=names)

    df["REOS"] = pass_list

    df = pd.concat([df, temp_df], axis=1)

    return df


def qed_filter(df):
    temp_df = []

    pass_list = []

    for item in df["smiles"]:
        molecule = Chem.MolFromSmiles(item)

        # Drug-Like (QED):
        #     mass < 400
        #     ring count > 0
        #     rotatable bond count < 5
        #     h-bond donor count <= 5
        #     h-bond acceptor count <= 10
        #     logP < 5

        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)
        h_bond_donor = Descriptors.NumHDonors(molecule)
        h_bond_acceptors = Descriptors.NumHAcceptors(molecule)
        rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
        num_of_rings = Chem.rdMolDescriptors.CalcNumRings(molecule)

        # print(molecular_weight,logp,h_bond_donor,h_bond_acceptors,formal_charge,rotatable_bonds,heavy_atoms)

        conditions = [
            molecular_weight <= 400,
            logp <= 5,
            h_bond_donor <= 5,
            h_bond_acceptors <= 10,
            num_of_rings > 0,
            rotatable_bonds <= 5,
        ]

        # descriptors = [num_of_rings]

        temp_df.append(num_of_rings)

        # print(conditions)
        pass_qed = conditions.count(True) == 6

        if pass_qed == True:
            pass_list.append(0)

        else:
            pass_list.append(1)

    df["QED"] = pass_list
    df["RingNum"] = temp_df

    return df


def calculate_filters(df, colname_chembl):
    df_smiles, unusedDrugs_df = GetSmiles(df, colname_chembl)

    # df_smiles = remove_salt(df_smiles,'smiles')

    temp = ro5_filter(df_smiles)
    temp = ghose_filter(temp)
    temp = veber_filter(temp)
    temp = reos_filter(temp)
    temp = qed_filter(temp)

    df = df.loc[df.groupby("drugId")["phase"].idxmax()]
    df = df.reset_index(drop=True)

    df = pd.merge(temp, df[["drugId", "phase"]], on="drugId", how="left")

    df[["phase"]] = df[["phase"]].astype(int)

    df = df.drop_duplicates()
    df = df.reset_index(drop=True)

    return (df, unusedDrugs_df)


def create_charts_with_calc_filters(df_filters_df):
    """
    This functions filters the dataframe obtained from calculate_filters(), and creates a stacked barchart and a piechart.
    """
    #    filter_cols = ["Lipinski_ro5", "Ghose", "Veber", "REOS", "QED"]
    #    df_filters = calc_filters_df[filter_cols]

    df_filters_df = df_filters_df.replace({0: "Yes", 1: "No"})

    grouped_stacked_filters = (
        df_filters_df.stack().groupby(level=[1]).value_counts().unstack()
    )

    stack_order = ["Yes", "No"]
    ordered_stacked_df = grouped_stacked_filters.reindex(
        columns=stack_order, fill_value=0
    )

    calc_fig = go.Figure()

    bar_colors = {
        "Yes": "#1f77b4",
        "No": "#ff7f0e",
    }  # Setting the colors exactly to match the github notebook color scheme

    for col in ordered_stacked_df.columns:
        calc_fig.add_trace(
            go.Bar(
                name=col,
                x=ordered_stacked_df.index,
                y=ordered_stacked_df[col],
                marker=dict(
                    color=bar_colors.get(col, "#CCCCCC"),
                    line=dict(width=2, color="white"),
                ),
            )
        )

    calc_fig.update_layout(
        barmode="stack",
        title="Bar plots showing drug-likeness across various filters",
        xaxis_title="Filters for drug-likeness",
        yaxis_title="Total number of drugs",
        width=800,
        height=700,
        plot_bgcolor="white",
        xaxis_tickangle=0,
        showlegend=True,
    )

    st.plotly_chart(calc_fig)
    return calc_fig


def create_pie_chart_from_calc_filters(calc_filters_df, df_filters):
    """
    This function creates a pie chart from the filtered dataframe obtained from calculate_filters().
    """
    #    filter_cols = ["Lipinski_ro5", "Ghose", "Veber", "REOS", "QED"]
    #    df_filters = calc_filters_df[filter_cols]

    calc_filters_df["Flag"] = df_filters.sum(axis=1, numeric_only=True)

    df_drugType = (
        calc_filters_df.groupby(["Flag"])["Flag"].count().reset_index(name="count")
    )

    df_drugType = df_drugType.rename(
        index={
            0: "Passed All Filters",
            1: "1 Violation",
            2: "2 Violations",
            3: "3 Violations",
            4: "4 Violations",
            5: "5 Violations",
        }
    )

    fig = go.Figure(
        data=[
            go.Pie(
                labels=df_drugType.index,
                values=df_drugType["count"],
                textinfo="label+percent",
                textposition="inside",
                marker=dict(
                    colors=[
                        "#1f77b4",
                        "#ff7f0e",
                        "#2ca02c",
                        "#d62728",
                        "#9467bd",
                        "#8c564b",
                    ],
                    line=dict(color="white", width=2),
                ),
            )
        ]
    )

    fig.update_layout(
        title="Pie Chart Showing the Number of Drugs Passing the Filters",
        width=800,
        height=700,
        plot_bgcolor="white",
    )

    st.plotly_chart(fig)
    return fig


def formulate_calc_filters_df(calc_filters_df, df_filters):
    """
    This function is used to restructure calc_filters dataframe before exporting to CSV.
    """
    #    filter_cols = ["Lipinski_ro5", "Ghose", "Veber", "REOS", "QED"]
    #    df_filters = calc_filters_df[filter_cols]
    calc_filters_df["Flag"] = df_filters.sum(axis=1, numeric_only=True)
    calc_filters_df.rename(
        columns={"smiles": "SMILES", "drugId": "DrugID"}, inplace=True
    )
    new_col_order = (
        ["DrugID", "SMILES"] + list(calc_filters_df.columns[2:-1]) + ["Flag"]
    )
    calc_filters_df = calc_filters_df[new_col_order]
    st.write(calc_filters_df)
    return calc_filters_df


def create_drug_likeness_zip(
    calc_figures, pie_figures, calc_filters_df, zip_filename="drug_likeness_results.zip"
):
    """
    Creates a zip file containing two interactive Plotly HTML charts and a XLSX file.
    """
    zip_buffer = io.BytesIO()

    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
        bar_html = calc_figures.to_html(full_html=False, include_plotlyjs="cdn")
        zip_file.writestr("drug_likeness_calc_figures.html", bar_html)

        pie_html = pie_figures.to_html(full_html=False, include_plotlyjs="cdn")
        zip_file.writestr("drug_likeness_pie_figures.html", pie_html)

        excel_buffer = io.BytesIO()
        with pd.ExcelWriter(excel_buffer, engine="openpyxl") as writer:
            calc_filters_df.to_excel(writer, index=False, sheet_name="Results")
        excel_buffer.seek(0)

        zip_file.writestr("drug_likeness_results.xlsx", excel_buffer.getvalue())

    zip_buffer.seek(0)
    return zip_buffer.getvalue()


def create_zip():
    drugs_df = state["drugs_df"].to_csv(index=False)
    dis2prot_df = state["dis2prot_df"].to_csv(index=False)
    dis2snp_df = state["dis2snp_df"].to_csv(index=False)
    advEff_df = state["adv_effect"].to_csv(index=False)
    humanProtDict = state["human_protein"]
    kg = state["graph"]
    if "viral_protein" in state:
        viralProtDict = state["viral_protein"]

    files = {
        "DiseaseAssociatedDrugs.csv": drugs_df,
        "DiseaseAssociatedProteins.csv": dis2prot_df,
        "DiseaseAssociatedSNPs.csv": dis2snp_df,
        "AdverseEffects.csv": advEff_df,
    }

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
        for file_name, file_content in files.items():
            zip_file.writestr(file_name, file_content)

        pickle_buffer_kg = io.BytesIO()
        pickle.dump(kg, pickle_buffer_kg)
        zip_file.writestr(f"{state.kg_name}.pkl", pickle_buffer_kg.getvalue())

        kg_csv_buffer = io.StringIO()
        pybel.to_csv(kg, kg_csv_buffer)
        zip_file.writestr(f"{state.kg_name}.csv", kg_csv_buffer.getvalue())

        kg_bel_buffer = io.StringIO()
        pybel.to_bel_script(kg, kg_bel_buffer)
        zip_file.writestr(f"{state.kg_name}.bel", kg_bel_buffer.getvalue())

        kg_graphml_buffer = io.BytesIO()
        pybel.to_graphml(kg, kg_graphml_buffer)
        zip_file.writestr(f"{state.kg_name}.graphml", kg_graphml_buffer.getvalue())

        pickle_buffer_humanProt = io.BytesIO()
        pickle.dump(humanProtDict, pickle_buffer_humanProt)
        zip_file.writestr(
            f"{state.kg_name}_UniProtDict.pkl",
            pickle_buffer_humanProt.getvalue(),
        )

        try:
            pickle_buffer_viralProt = io.BytesIO()
            pickle.dump(viralProtDict, pickle_buffer_viralProt)
            zip_file.writestr(
                f"{state.kg_name}_UniProtDict.pkl",
                pickle_buffer_viralProt.getvalue(),
            )
        except NameError:
            viralProtDict = None

        if "figures" in state:
            for fig_name, fig in state["figures"].items():
                fig_buffer = io.BytesIO()
                fig.write_image(fig_buffer, format="png", engine="kaleido")
                fig_buffer.seek(0)
                zip_file.writestr(f"{fig_name}.png", fig_buffer.getvalue())

    # zip_buffer.seek(0)

    return zip_buffer.getvalue()


def load_pickle_file(uploaded_file):
    """
    Loads the uploaded .pkl file and returns the BELGraph object if valid.
    """
    try:
        # uploaded_file is already a file-like object
        query_graph = pickle.load(uploaded_file)

        # Check if it's a valid BELGraph
        if not isinstance(query_graph, pybel.BELGraph):
            st.error("Uploaded pickle file does not contain a BELGraph object.")
            st.stop()

        return query_graph

    except Exception as e:
        st.error(f"Failed to load pickle file: {e}")
        st.stop()


def query_graph_info(graph_data):
    """
    Reads the graph and displays information about the graph.
    """
    total_nodes = graph_data.number_of_nodes()
    total_edges = graph_data.number_of_edges()
    total_relations = Counter(
        data["relation"] for _, _, data in graph_data.edges(data=True)
    )

    # Create a manual summary if summarize() doesn't work
    summary = (
        f"Nodes: {total_nodes}, Edges: {total_edges}, Relations: {len(total_relations)}"
    )

    st.markdown(f"### Summary of your graph:\n{summary}")
    st.markdown(f"**Total no. of relations: {len(total_relations)}")
    st.markdown(f"**Total no. of edges: {total_edges}")
    st.markdown(f"**Total no. of nodes: {total_nodes}")
    st.markdown("**Relation types:**")
    for rel, count in total_relations.most_common():
        st.markdown(f"- {rel}: {count}")


def display_interactive_belgraph(graph_data):
    graph_subset = induction.get_random_subgraph(graph_data)

    total_nodes = graph_subset.number_of_nodes()
    total_edges = graph_subset.number_of_edges()
    st.markdown(
        f"Taking a random subgraph with **{total_nodes} nodes** and **{total_edges} edges...**"
    )
    graph_subset_html = pybel_jupyter.to_html(graph_subset)
    st.components.v1.html(graph_subset_html, height=600, scrolling=True)

    return graph_subset_html


def download_interactive_belgraph(graph_html):
    """Lets the user download the BELGraph visualization as an HTML file."""

    if not isinstance(graph_html, str):
        st.error("Error: The graph visualization must be generated first.")
        return

    html_bytes = graph_html.encode("utf-8")
    buffer = io.BytesIO(html_bytes)

    st.download_button(
        label="Download BELGraph as HTML",
        data=buffer,
        file_name="belgraph_visualization.html",
        mime="text/html",
        help="Download the HTML format of this interactive graph.",
        on_click="ignore",
    )
