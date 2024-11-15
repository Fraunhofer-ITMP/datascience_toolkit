"""App for Knowledge Graph Generator (KGG)"""

import base64
import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(
    layout="wide",
    page_title="PROTACKB",
    page_icon=":books:",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
        button[title^=Exit]+div 
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    """
        <style>
            .block-container {
                padding-top: 1.5rem;
                padding-bottom: 1.5rem;
                padding-left: 5rem;
                padding-right: 5rem;
            }
            .stTabs [data-baseweb="tab-list"] button [data-testid="stMarkdownContainer"] p {font-size:1.3rem;}
        </style>
        """,
    unsafe_allow_html=True,
)  # .block-conatiner controls the padding of the page, .stTabs controls the font size of the text in the tabs

st.markdown(
    "<h1 style='text-align: center; color: #149372;'> PROTACKB</h1> <br>",
    unsafe_allow_html=True,
)

st.header("Description", divider="gray")
st.write(
    """**PROTACKB** is a comprehensive knowledgebase of PROTACs, molecular glues and degraders developed with the motive to harmonize and consolidate relevant data from public databases, patents and ProxiDrugs project. As of now, it consists of 14 types of entities (**24K unique entities**) and **55K relationships** which are represented using graph theory based data structure model. The entities are enriched with knowledge from curated databases such as Ubinet, OpenTargets, UniProt and ChEMBL. The PROTACKB is hosted in Neo4j, a visualization tool with intuitive graphical user interface (GUI). The access to Neo4j server is through personalized login credentials and is only available to members of ProxiDrugs consortium. Additional information about ProxiDrugs can be found [here](https://www.proxidrugs.de/)."""
)

video_file = open("./images/pages/ARdata_protackb.mp4", "rb")
video_bytes = video_file.read()

st.video(video_bytes, loop=True, start_time="8s", end_time="48s", autoplay=True)

st.header("General workflow", divider="gray")

kg_worklow = "./images/pages/workflow_ptackb.png"
st.image(
    kg_worklow,
    caption="A schematic representation of PROTACKB",
    width=1000,
)

st.write(
    """The process of creating PROTACKB involves various steps for normalizing and standardizing chemical representations of entities such as PROTACs, linkers, warheads, e3 binders, glues and degraders. These include conversion of SDF structures to SMILES/InChI, removal of salts from structures, removal of invalid structures, generation of canonical SMILES and InChI key and finally assignment of unique identifiers (CBDREGNo.) to entities. Afterwards, using SMILES/InChI Key as a query to ChEMBL and OpenTargets API, the knowledge about mechanism of action, associated disease and adverse effect is retrieved (if any). Similary, proteins are first mapped to UniProt database and enriched with the knowledge of associated pathway, molecular function and biological process using UniProt's API. The e3 ligases are enriched with knowledge from UbiNet which are about classes of e3 ligases and domain structure. This knowledge is represented in the form of triples (i.e., :blue[subject]-:red[_predicate_] -:green[object] which are the building blocks of PROTACKB. For Examples: :blue[**PROTAC123**] :red[**_hasWarhead_**] :green[**Warhead456**] and :blue[**Warhead456**] :red[**_Targets_**] :green[**TYK2**]. The entity specific metadata is added to the corresponding entity, for instance: chemical structures will always have descriptors such as Mol Wt, TPSA, HBA, HBD, etc. and protein will mostly have UniProt id and links to external databases. The PROTACKB entities are also cross-linked to the PROXIDRUGSDB which is a database for assay and experimental data generated in ProxiDrugs."""
)

st.header("Overview of entities in PROTACKB", divider="gray")

path_ns = "./images/pages/KG_namespace.html"
HtmlFile = open(path_ns, "r", encoding="utf-8")
source_code = HtmlFile.read()
components.html(source_code, height=450)

col1, col2 = st.columns(2)

with col1:
    st.write(
        """
        **What is E3 ligase?**: A ubiquitin ligase (also called an E3 ubiquitin ligase) is a protein that recruits an E2 ubiquitin-conjugating enzyme that has been loaded with ubiquitin, recognizes a protein substrate, and assists or directly catalyzes the transfer of ubiquitin from the E2 to the protein substrate. (Ref: [Wikipedia](https://en.wikipedia.org/wiki/Ubiquitin_ligase))

        **What do we know?**: There are about 400 known human E3 ligases, however only a handful of E3 ligases are currently used in PROTAC-design. Following figure shows the proportion of 1224 E3 binders for each E3 ligase."""
    )

with col2:
    path_ns = "./images/pages/E3Ligase.html"
    HtmlFile = open(path_ns, "r", encoding="utf-8")
    source_code = HtmlFile.read()
    components.html(source_code, height=450)


st.header("What can you do with PROTACKB?", divider="gray")

st.markdown(
    """<div style="text-align: justify;"> The main purpose of PROTACKB is to systematically assemble information and knowledge about Targeted Protein Degradation (TPD), especially PROTACs and their building blocks. In this regard, PROTACKB is a megastore as it provides information about 14 different types of entities, all of which are related to TPDs. The Neo4j server, which hosts the PROTACKB, provides a search interface to look for entities of your interest and know more about it. For instance, if you are interested in a PROTAC, you can get information about its physicochemical properties, standard representations such as SMILES and InChI key, its provenance (link to scientific article, PubMed or patent number) and so on. Well, these are just the metadata about the PROTAC. You can eventually expand the neighborhood around it to identify its correpsonding linker, warhead, e3 binder and target protein. Once done, you can expand the neighborhood around the target protein to know which other PROTACs are also out there or get to know the e3 ligase for e3 binder and PROTAC. These are simple things you can easily do with just a couple of clicks. Start with a entity of interest and get to know more. And all of this can be saved for your next visit to Neo4j, you don't have to restart. There are features also to export screenshots and CSV files of your current instance. </div>""",
    unsafe_allow_html=True,
)

st.markdown("<br>", unsafe_allow_html=True)

st.markdown(
    """<div style="text-align: justify;"> There's more! Going after entity of your interest is straightforward and might only reveal just a little more than what you might already know. Can we possibly ask the PROTACKB to find something more complex in the mix of 24K nodes and 55K relationships? Absolutely! That's why Neo4j comes with Cypher. It is a query language to formulate complex questions with all your lists of this and that but not that. For instance, you can write a query to know if warheads are used as e3 binders or vice-versa in PROTACs or identify PROTAC(s) with identical e3 binder and warhead. That would be a lazy design if it exists! Something more complicated but realistic would be to identify PROTAC(s) and the target protein for a certain disease. And because you know that the affect organ expresses a certain e3 ligase abundantly, you can easily add this condition to retrieve PROTAC(s)/e3 binder(s) with only this e3 ligase. Lastly, you would also like the PROTAC(s) to be registered in ChEMBL. Not a problem at all except one. Since Cypher is a language of its own, one needs to formulate this question as a Cypher query. That is difficult but worry not, ChatGPT and Perplexity like AI tools can easily do this. In the near future, a ChatGPT integrated Neo4j will formulate your natural language query to a Cypher query. Just about time, we are working on it. </div>""",
    unsafe_allow_html=True,
)

st.markdown("<br>", unsafe_allow_html=True)

st.write(
    """<div style="text-align: justify;"> Lastly, not everything has to be about PROTACKB. It is just a comprehensive resource with TPD knowledge. The source data about TPD counterparts can be used in many different ways such as developing ML/AI tools for characterization, profiling and virtual screening. Define your use case and applications and just write to us. </div>""",
    unsafe_allow_html=True,
)
st.markdown("<br>", unsafe_allow_html=True)

# footer with text and green background
st.markdown(
    "<footer style='background-color: #149372; padding: 10px; border-radius: 10px;'>"
    "<p style='color: white; text-align: center;'>Fraunhofer ITMP © 2024</p>"
    "<p style='color: white; text-align: center;'>This work has been conducted across several key projects in which ITMP has been actively involved.</p>"
    "</footer>",
    unsafe_allow_html=True,
)
