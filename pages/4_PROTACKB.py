import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(
    layout="wide",
    page_title="PROTACKB",
    page_icon="🕸️",
    initial_sidebar_state="collapsed",
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

st.markdown("### Description")
st.write("**PROTACKB** is a comprehensive knowledgebase of PROTACs, molecular glues and degraders developed with the motive to harmonize and "
         "consolidate relevant data from public databases, patents and ProxiDrugs project. As of now, it consists of 15 types of entities "
         "(**27K unique nodes**) and **115K relationships** which are represented using graph theory based data structure model. The entities are enriched with knowledge from curated databases such as Ubinet, OpenTargets, UniProt "
         "and ChEMBL. The PROTACKB is hosted in Neo4j, a visualization tool with intuitive graphical user interface (GUI). The access to Neo4j server is through "
         "personalized login credentials and is only available to members of ProxiDrugs consortium. Additional information about ProxiDrugs can be found [here](https://www.proxidrugs.de/).")

st.markdown("### General workflow")
st.write("The process of creating PROTACKB involves various steps for normalizing and standardizing chemical representations of entities such as PROTACs, linkers, warheads, "
         "e3 binders, glues and degraders. These include conversion of SDF structures to SMILES/InChI, removal of salts from structures, removal of invalid structures, "
         "generation of canonical SMILES and InChI key and finally assignment of unique identifiers (CBDREGNo.) to entities. Afterwards, using SMILES/InChI Key as a query to "
         "ChEMBL and OpenTargets API, the knowledge about mechanism of action, associated disease and adverse effect is retrieved (if any). Similary, proteins are first mapped to "
         "UniProt database and enriched with the knowledge of associated pathway, molecular function and biological process using UniProt's API. The e3 ligases are enriched with "
         "knowledge from UbiNet which are about classes of e3 ligases and domain structure. This knowledge is represented in the form of triples (i.e., :blue[subject]-:red[_predicate_]"
         "-:green[object] which are the building blocks of PROTACKB. For Examples: :blue[**PROTAC123**] :red[**_hasWarhead_**] :green[**Warhead456**] and :blue[**Warhead456**] :red[**_Targets_**] :green[**TYK2**]. The entity specific metadata is added to the corresponding entity, for instance: "
         "chemical structures will always have descriptors such as Mol Wt, TPSA, HBA, HBD, etc. and protein will mostly have UniProt id and links to external databases. The PROTACKB entities "
         "are also cross-linked to the PROXIDRUGSDB which is a database for assay and experimental data generated in ProxiDrugs.")

kg_worklow = "./images/workflow_ptackb.png"
st.image(
    kg_worklow,
    caption="A schematic representation of PROTACKB",
)

st.markdown("### PROTACKB Demo: Visualization of Androgen Receptor (AR) protein in Neo4j")

video_file = open("./images/ARdata_protackb.mp4","rb")
video_bytes = video_file.read()

st.video(video_bytes,loop=True,start_time="8s",end_time="48s",autoplay=True)

st.markdown("### PROTACKB statistics")

st.write("In the following figures and plots, we provide various insights of the PROTACKB.")

st.markdown("##### Overview of entities in PROTACKB")

path_ns = './images/KG_namespace.html'
HtmlFile = open(path_ns, 'r', encoding='utf-8')
source_code = HtmlFile.read()
# hit and trial for height adjustment == 400
components.html(source_code, height=450, width=800)

st.markdown("##### E3 ligases and relative percentage of E3 binders")

st.write('There are about 400 known human E3 ligases, however only a handful of E3 ligases are currently used in PROTAC-design. '
         'Following figure shows the proportion of 1224 E3 binders for each E3 ligase.')

path_ns = './images/E3Ligase.html'
HtmlFile = open(path_ns, 'r', encoding='utf-8')
source_code = HtmlFile.read()
# hit and trial for height adjustment == 400
components.html(source_code, height=450, width=800)




