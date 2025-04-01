# -*- coding: utf-8 -*-

import os
import pickle
import logging
import json
import pandas as pd
import networkx as nx
import requests
from collections import defaultdict
import pybel

import zipfile
import io

import streamlit as st
from stqdm import stqdm
from pybel import BELGraph
from pybel.dsl import Protein, Abundance, BiologicalProcess, Pathology, Gene

from pandasgwas import get_variants
from pandasgwas.get_variants import get_variants_by_efo_id

from chembl_webresource_client.new_client import new_client

from rdkit import Chem
from rdkit.Chem import Descriptors

import plotly.graph_objects as go
from plotly.subplots import make_subplots

pd.options.mode.chained_assignment = None  # default='warn'


logger = logging.getLogger("__name__")

DATA_DIR = "data/"


def load_kg(path):
    infile = open(path, "rb")
    kg = pickle.load(infile)
    infile.close()
    return kg


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
    if "graph_stats_plot" in st.session_state:
        st.session_state["graph_stats_plot"] = fig

    st.session_state["graph_stats_plot"] = fig

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

    if "graph_ns_plot" in st.session_state:
        st.session_state["graph_ns_plot"] = figure_ns
    st.session_state["graph_ns_plot"] = figure_ns

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


def GetDiseaseAssociatedDrugs(disease_id, CT_phase) -> pd.DataFrame:
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


def GetDiseaseAssociatedProteins(
    disease_id, index_counter=0, merged_list=[]
) -> pd.DataFrame:

    efo_id = str(disease_id)

    query_string = """
        query associatedTargets($efoId: String!,$index:Int!){
          disease(efoId: $efoId){
            id
            name
            associatedTargets(page:{size:3000,index:$index}){
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
    # query_string = query_string.replace("$efoId",f'"{efo_id}"')

    variables = {"efoId": efo_id, "index": index_counter}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    # r = requests.post(base_url, json={"query": query_string})
    # print(r.status_code)

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    result = api_response["data"]["disease"]["associatedTargets"]["rows"]

    merged_list.extend(result)

    if result:
        counter = index_counter + 1
        GetDiseaseAssociatedProteins(disease_id, counter, merged_list)

    temp_list = []
    for item in merged_list:
        # api_response['data']['disease']['associatedTargets']['rows']
        # print(item['target'])
        # break
        for obj in item["target"]["proteinIds"]:
            if obj["source"] == "uniprot_swissprot":
                # print(obj)
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

    df.drop_duplicates(inplace=True)  # recursion function created duplicates
    return df


def GetDiseaseAssociatedProteinsPlot(df) -> None:
    """Plotting the protein confidence scores associated with a disease."""
    st.markdown("**Protein-Disease Association summary**")
    disease_id = df["disease_id"].unique()[0]

    assert df["disease_id"].nunique() == 1

    st.markdown(
        f"We have identified :green[{len(df)}] proteins (Swiss-Prot) associated with the disease (:green[{disease_id}]). "
        + "Following is a histogram that shows distribution of proteins based on scores provided by OpenTargets. "
        + "The scores are influenced by various factors such as genetic associations, expression, mutations, known pathways, "
        + "targeting drugs and so on."
    )

    df.sort_values(by="Score", ascending=False, inplace=True)

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
    st.markdown(
        "> **_NOTE:_** 📝 Please note that the proteins identified may not be unique if you combined two or more diseases."
    )


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


def GetDiseaseSNPs(disease_id) -> pd.DataFrame:
    try:
        snps = get_variants_by_efo_id(disease_id)
        snps_df = snps.genomic_contexts
        snps_df["disease_id"] = disease_id
        snps_df = snps_df.reset_index(drop=True)
    except:
        st.markdown(
            "> <i>No SNPs found. This could be either because: (1.) The identifier of your disease of interest is not compatible with Experimental Factor Ontology (EFO) or (2.) No SNPs have been reported for the disease.</i>",
            unsafe_allow_html=True,
        )
        snps_df = pd.DataFrame()

    return snps_df


def GetViralProteins(query_disease) -> list:
    # file downloaded from https://www.genome.jp/ftp/db/virushostdb Dated: 12/09/2023
    virus = pd.read_csv(
        "https://raw.githubusercontent.com/Fraunhofer-ITMP/kgg/main/data/misc/virushostdb.csv"
    )

    cols = ["virus tax id", "virus name", "DISEASE", "host tax id"]
    virus = virus[cols]

    virus = virus.rename(
        columns={
            "virus tax id": "virus tax id",
            "virus name": "virus name",
            "DISEASE": "disease",
            "host tax id": "host tax id",
        }
    )

    # filter virus with host humans
    virus = virus.loc[virus["host tax id"] == 9606.0, :]
    virus = virus.reset_index(drop=True)

    # replace 9606.0 to 9606
    virus["host tax id"] = pd.to_numeric(virus["host tax id"], downcast="integer")

    # subset df with disease keyword
    virus_subset_1 = virus[
        virus["disease"].str.contains(query_disease, na=False, case=False)
    ]

    if virus_subset_1.empty:
        return []

    st.warning(
        "The workflow has identified your query as a viral disease. In the next steps, you can choose to add viral proteins (SWISS-Prot) in the KG. \n"
    )

    # Cleaning the dataframe for outputting
    virus_subset_1.drop(columns=["host tax id"], inplace=True)
    virus_subset_1.drop_duplicates(inplace=True)

    viral_output_req = st.text_input(
        "Do you want to add viral proteins to the KG?",
        value="No",
    )

    if viral_output_req.lower() == "no":
        st.write("You can skip this step and proceed to the next step.")
        uprot_list = []

    elif viral_output_req.lower() == "yes":

        st.markdown(
            "**The following viruses have been identified affecting humans, please search and find the virus of interest.**"
        )
        config = {
            "index": st.column_config.NumberColumn("Index"),
            "virus name": st.column_config.TextColumn("Virus Name"),
            "disease": st.column_config.TextColumn("Disease"),
            "host tax id": st.column_config.TextColumn("Host Tax ID"),
        }
        st.dataframe(
            virus_subset_1,
            use_container_width=True,
            column_config=config,
        )

        first_index_of_subset = virus_subset_1.index[0]
        temp_id = st.text_input(
            "Enter the index value(s). If multiple, use space, for example -> 0 1 3: ",
            value=str(first_index_of_subset),
        )

        temp_ids = temp_id.split(" ")
        temp_id_list = [int(x) for x in temp_ids]

        uprot_list = []

        for item in temp_id_list:
            tax_id = virus_subset_1.loc[item]["virus tax id"]

            # fetch tax id related proteins from Uniprot
            # the link can be created from downloads option in uniprot
            query_string = (
                "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Clength%2Cgene_primary%2Cprotein_name&format=tsv&query=%28%28taxonomy_id%3A"
                + str(tax_id)
                + "%29+AND+%28reviewed%3Atrue%29%29"
            )

            query_uniprot = requests.get(query_string)
            query_uniprot = query_uniprot.text.split("\n")

            query_uniprot_df = pd.DataFrame(
                [x.strip().split("\t") for x in query_uniprot]
            )
            cols = query_uniprot_df.iloc[0]

            query_uniprot_df = query_uniprot_df[1 : len(query_uniprot_df) - 1]
            query_uniprot_df.columns = cols
            temp = list(query_uniprot_df["Entry"])

            uprot_list.append(temp)

        uprot_list = [item for sublist in uprot_list for item in sublist]

    st.write(
        f"A total of :green[{str(len(uprot_list))} viral proteins] have been identified."
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


def createInitialKG(efo_id, ct_phase):
    """Creating the initial Knowledge Graph using the disease and protein data."""
    for functions in stqdm(
        ["disease_drugs", "disease_proteins", "disease_snp"],
        "Fetching real-time data from databases. Be patient!",
    ):
        if functions == "disease_drugs":
            drugs_df = GetDiseaseAssociatedDrugs(efo_id, ct_phase)
            if not drugs_df.empty:
                drugs_df = drugs_df.reset_index(drop=True)

        elif functions == "disease_proteins":
            dis2prot_df = GetDiseaseAssociatedProteins(efo_id)
            if not dis2prot_df.empty:
                dis2prot_df = dis2prot_df.reset_index(drop=True)

        elif functions == "disease_snp":
            dis2snp = GetDiseaseSNPs(efo_id)
            if not dis2snp.empty:
                dis2snp = dis2snp.reset_index(drop=True)

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


def GetAdverseEvents(chem_list) -> pd.DataFrame:
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
                ): "https://www.ebi.ac.uk/chembl/compound_report_card/"
                + item
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
                ): "https://www.ebi.ac.uk/QuickGO/term/"
                + gobp_id
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
                ): "https://www.ebi.ac.uk/QuickGO/term/"
                + gomf_id
            },
            "QuickGO",
        )
        nx.set_node_attributes(
            graph,
            {BiologicalProcess(namespace="GOMF", name=item): gomf_id},
            "Gene Ontology identifier",
        )

    return graph


def finalizeKG(
    filtered_protein_df: pd.DataFrame,
    drugs_df: pd.DataFrame,
    kg_name: str,
    dis2snp_df: pd.DataFrame,
    viral_proteins: list,
):
    """Finalizing the Knowledge Graph by adding proteins and drugs data."""

    # create empty KG
    kg = BELGraph(name=kg_name, version="0.0.1")

    # Prots
    unique_proteins = list(set(filtered_protein_df["UniProt"]))

    uprot_ext = ExtractFromUniProt(unique_proteins)
    kg = uniprot_rel(uprot_ext, "HGNC", kg)
    kg = gene_ontology_annotation(kg, uprot_ext)

    if len(viral_proteins) > 0:
        vir_uprot_ext = ExtractFromUniProt(viral_proteins)
        kg = uniprot_rel(vir_uprot_ext, "VP", kg)

    # Drugs
    if not drugs_df.empty:
        st.write(
            "A total of "
            + str(len(list(set(drugs_df["drugId"]))))
            + " drugs have been identified. Now fetching relevant data"
        )

        chembl2mech = RetMech(list(set(drugs_df["drugId"])))
        chembl2act = RetAct(list(set(drugs_df["drugId"])))

        prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(
            chembl2mech
        )
        prtn_as_chembl = set(prtn_as_chembl)
        prtn_as_chembl = list(prtn_as_chembl)
        chembl2uprot = chembl2uniprot(prtn_as_chembl)

        chembl2act_2 = chembl2gene2path(chembl2uprot, chembl2act)
        chembl2mech_2 = chembl2gene2path(chembl2uprot, chembl2mech)

        kg = chem2moa_rel(chembl2mech_2, "HGNC", unique_proteins, kg)
        kg = chem2act_rel(chembl2act_2, "HGNC", unique_proteins, kg)
        kg = gene2path_rel(chembl2uprot, "HGNC", unique_proteins, kg)

        adv_effect = GetAdverseEvents(list(set(drugs_df["drugId"])))

        kg = chembl2adverseEffect_rel(adv_effect, kg)

        kg = chembl_name_annotation(kg, drugs_df)
    else:
        st.write("No drugs have been identified for the disease.")
        adv_effect = pd.DataFrame()

    # SNPs
    if not dis2snp_df.empty:
        kg = snp2gene_rel(dis2snp_df, kg)

    st.write("Your KG is now generated!", "\n")
    return kg, adv_effect


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


def refactor_state(session_keys: list):
    # Delete the keys from session state
    for key in session_keys:
        if key in st.session_state:
            del st.session_state[key]


def save_files():
    folder_name = st.session_state.disease_name
    os.makedirs(folder_name, exist_ok=True)

    # Save DataFrames as CSV files
    files = {
        "DiseaseAssociatedDrugs.csv": st.session_state["drugs_df"],
        "DiseaseAssociatedProteins.csv": st.session_state["dis2prot_df"],
        "DiseaseAssociatedSNPs.csv": st.session_state["dis2snp_df"],
        "AdverseEffects.csv": st.session_state["adv_effect"],
    }

    for file_name, df in files.items():
        df.to_csv(os.path.join(folder_name, file_name), index=False)

    # Save other data structures
    kg = st.session_state["graph"]
    # humanProtDict = st.session_state["human_protein"]
    viralProtDict = st.session_state["viral_prot"]

    bel_file = os.path.join(folder_name, f"{st.session_state.kg_name}.bel")
    graphml_file = os.path.join(folder_name, f"{st.session_state.kg_name}.graphml")

    pybel.dump(kg, bel_file)
    pybel.to_graphml(kg, graphml_file)

    if viralProtDict is not None:
        viral_pickle_file = os.path.join(
            folder_name, f"{st.session_state.kg_name}_ViralUniProtDict.pkl"
        )
        pickle.dump(viralProtDict, open(viral_pickle_file, "wb"))

    # Save the plotly figures
    ns_fig = st.session_state["graph_ns_plot"]
    stat_fig = st.session_state["graph_stats_plot"]
    ns_fig.write_image(
        os.path.join(folder_name, f"{st.session_state.kg_name}_NamespaceSummary.png")
    )
    stat_fig.write_image(
        os.path.join(folder_name, f"{st.session_state.kg_name}_GraphStatistics.png")
    )

    # Create zip file
    zip_path = f"{folder_name}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(folder_name):
            for file in files:
                file_path = os.path.join(root, file)
                zipf.write(file_path, os.path.relpath(file_path, folder_name))

    return zip_path
