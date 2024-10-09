"""App for Knowledge Graph Generator (KGG)"""

import requests
import streamlit as st
import streamlit.components.v1 as components

from utils import disease_figures


st.set_page_config(
    layout="wide",
    page_title="Knowledge Graph Generator (KGG)",
    page_icon="🕸️",
    initial_sidebar_state="collapsed",
)

st.markdown(
    "<h1 style='text-align: center; color: #149372;'>Knowledge Graph Generator (KGG)</h1> <br>",
    unsafe_allow_html=True,
)

tab1, tab2 = st.tabs(
    [
        "Description",
        "DisKGG",
    ]
)


with tab1:
    st.markdown("### What is KGG?")
    st.markdown(
        """**K**nowledge **G**raph **G**enerator (KGG) is a fully automated workflow for creating disease-specific KGs. It is developed for a broad spectrum of researchers and scientists, especially for those who are into pre-clinical drug discovery, understanding disease mechanisms/comorbidity and drug-repurposing. The KGG embeds underlying schema of curated public databases to retrieve relevant knowledge which is regarded as the gold standard for high quality data. The KGG is leveraged on our previous contributions to the BY-COVID project where we developed workflows for identification of bio-active analogs for fragments identified in COVID-NMR studies (Berg, H et al.) and representation of Mpox biology (Karki, R et al.). The programmatic scripts and methods for KGG are written in python (version 3.10) and are available on [GitHub](https://github.com/Fraunhofer-ITMP/kgg).
        """
    )

    st.markdown(
        "Additional information about KGG can be found [here](https://www.infectious-diseases-toolkit.org/showcase/knowledge-graph-generator)"
    )

    kg_worklow = "./images/KGG_workflow_updated.png"
    st.image(
        kg_worklow,
        caption="A schematic representation of the 3 phases of KGG workflow.",
    )

    st.markdown("### Examples and Usecases")

    option = st.selectbox(
        label="Please choose a disease by clicking on the grey area below",
        options=("Alzheimer", "Parkinson", "Depression", "COVID-19"),
        index=0,
    )

    disease_figures(option)

with tab2:

    def getImage():

        img = "https://github.com/Fraunhofer-ITMP/kgg/blob/main/data/KGG.png?raw=true"
        response = requests.get(img, stream=True)
        img = Image.open(response.raw)

        return img


def getQuery():

    st.image(getImage())

    st.markdown(
        "**:blue[Welcome to the KG Generator tool. In the following steps, we will need some inputs from your side.]**"
    )

    with st.form(key="search_disease"):
        query = st.text_input(
            "Please enter the disease you are interested in and we will try to find the best matches for you."
        )

        submit_button = st.form_submit_button(label="Perform search")

    # st.write(st.session_state)

    if submit_button:

        st.success("Query received, now searching...", icon="✅")
        st.session_state.query = query

        return query

    return None

    # st.write(st.session_state)


def disease_df(query):

    # if 'query' not in st.session_state:
    #     st.session_state.query = None

    if query:

        search_df = searchDisease(query)

        if not search_df.empty:
            st.write(
                "**:blue[Here you go! Hopefully your disease of interest is in the list. If so, let's get started.]**"
            )

            st.session_state.disease = query
            st.session_state.search_df = search_df
            st.session_state.query_generated = True
            return True

        else:
            st.write("Ooops!! Did you have a typo in the name. Please try again!")
            st.session_state.query_generated = False

            return False


def get_inputs():

    # if 'efo_id' not in st.session_state:
    #     st.session_state.efo_id = None
    # if 'ctphase' not in st.session_state:
    #     st.session_state.ctphase = None
    # if 'kg_name' not in st.session_state:
    #     st.session_state.kg_name = None
    # if 'query_phaseCompleted' not in st.session_state:
    #     st.session_state.query_phaseCompleted = False

    with st.form(key="get_disease"):

        efo_id = st.text_input(
            "Please enter the index value of your disease of interest.", value=0
        )
        submit_disease = st.form_submit_button(label="Submit")

    if submit_disease:
        st.success("Index saved", icon="✅")
        st.session_state.efo_id = efo_id
        # st.write(st.session_state.efo_id)

    with st.form(key="choose_ct_phase"):

        st.write(
            "Please enter the clinical trial phase of chemicals which should be identified by the workflow. Use a number between 1 (early phase) and 4 (FDA approved). For example, if you use 3, the KG will fetch chemicals that are in phase 3. Also, remember that lower the input value, higher will be the number of identified chemicals and therefore the running time of workflow also increases."
        )

        ctphase = st.text_input("Your desired clinical trial phase:", value=4)
        submit_ctphase = st.form_submit_button(label="Submit")

    if submit_ctphase:
        st.success("Clinical trial phase saved", icon="✅")
        # st.write(submit_ctphase)
        st.session_state.ctphase = ctphase

    with st.form(key="nameKG"):

        kg_name = st.text_input("Please provide a name for you KG.", value="test_kg")
        submit_kgName = st.form_submit_button(label="Submit")

    if submit_kgName:
        st.success("KG name saved", icon="✅")
        # st.write(submit_kgName)
        st.session_state.kg_name = kg_name

    if (
        st.session_state.efo_id is not None
        and st.session_state.ctphase is not None
        and st.session_state.kg_name is not None
    ):
        st.success("All 3 inputs saved", icon="✅")
        st.session_state.query_phaseCompleted = True
        return (
            st.session_state.efo_id,
            st.session_state.ctphase,
            st.session_state.kg_name,
        )

    else:
        st.session_state.query_phaseCompleted = False
        return None


# if 'virus_search_submitted' not in st.session_state:
#     st.session_state.virus_search_submitted = False

# def getViralProteins(query_disease):
#     # file downloaded from https://www.genome.jp/ftp/db/virushostdb Dated: 12/09/2023
#     virus = pd.read_csv('https://raw.githubusercontent.com/Fraunhofer-ITMP/kgg/main/data/virushostdb.csv')
#
#     cols = ['virus tax id', 'virus name', 'DISEASE', 'host tax id']
#     virus = virus[cols]
#
#     # st.write(virus)
#
#     # filter virus with host humans
#     virus = virus.loc[virus['host tax id'] == 9606.0, :]
#     virus = virus.reset_index(drop=True)
#
#     # replace 9606.0 to 9606
#     virus["host tax id"] = pd.to_numeric(virus["host tax id"], downcast='integer')
#
#     # subset df with disease keyword
#     virus_subset_1 = virus[virus['DISEASE'].str.contains(query_disease, na=False, case=False)]
#
#     if not virus_subset_1.empty:
#
#         st.write('not empty')
#         st.write('\n')
#         st.write(
#             'The workflow has identified your query as a viral disease. Its proteins (SWISS-Prot) will be now represented in the KG.',
#             '\n')
#
#         st.write('Enter Viral Prot Function')
#
#         time.sleep(0.1)
#
#         # Initialize session state for virus name if not exists
#         if 'virus_name' not in st.session_state:
#             st.session_state.virus_name = ''
#         if 'virus_search_submitted' not in st.session_state:
#             st.session_state.virus_search_submitted = False
#
#         with st.form(key='identify_viral_proteins'):
#
#             virus_name = st.text_input("Do you want to look further for a specific virus? Please type its name or skip it by typing \'no\': ",
#                                        value=st.session_state.virus_name)
#
#             submit_button = st.form_submit_button(label='Search')
#
#         if submit_button:
#             #st.write(st.session_state.virus_name, ' submit button')
#             st.session_state.virus_name = virus_name
#             st.session_state.virus_search_submitted = True
#
#             if st.session_state.virus_search_submitted:
#                 st.write(st.session_state.virus_name, ' submit button')
#
#             # subset df with virus name
#             if st.session_state.virus_name.lower() != 'no':
#
#                 virus_subset_2 = virus[virus['virus name'].str.contains(virus_name, na=False, case=False)]
#                 #st.write(virus_subset_2)
#                 if not virus_subset_2.empty:
#                     st.write("Virus subset based on your input:")
#                     st.write(virus_subset_2)
#                 else:
#                     st.write("No matching viruses found for the given name.")
#                 # break
#             else:
#                 virus_subset_2 = pd.DataFrame()
#
#                 # merge subsets of df_1 and df_2
#                 virus_subset_merge = pd.concat([virus_subset_1, virus_subset_2])
#
#                 virus_subset_merge = virus_subset_merge.drop_duplicates(keep='first')
#
#                 virus_subset_merge = virus_subset_merge.reset_index(drop=True)
#
#                 virus_subset_merge['index'] = virus_subset_merge.index
#
#                 virus_subset_merge.style.hide(axis='index')
#
#                 virus_subset_merge = virus_subset_merge[['index', 'virus tax id', 'virus name', 'DISEASE', 'host tax id']]
#                 st.write(virus_subset_merge)
#
#                 #display(HTML(virus_subset_merge.to_html(index=False)))
#
#                 time.sleep(0.1)
#
#                 temp_id = st.text_input('Enter the index value(s). If multiple, use space, for example -> 0 1 3: ')
#
#                 #temp_id = input('Enter the index value(s). If multiple, use space, for example -> 0 1 3: ')
#
#                 st.write('\n')
#
#                 temp_id = temp_id.split(' ')
#                 temp_id = [int(x) for x in temp_id]
#                 # st.write(virus_subset_merge.loc[0]['virus tax id'])
#
#                 uprot_list = []
#
#                 for item in temp_id:
#                     tax_id = virus_subset_merge.loc[item]['virus tax id']
#                     # st.write(tax_id)
#
#                     # fetch tax id related proteins from Uniprot
#                     # the link can be created from downloads option in uniprot
#                     query_string = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Clength%2Cgene_primary%2Cprotein_name&format=tsv&query=%28%28taxonomy_id%3A' + str(
#                         tax_id) + '%29+AND+%28reviewed%3Atrue%29%29'
#
#                     # query_string = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Clength%2Cgene_primary%2Cprotein_name&format=tsv&query=%28%28taxonomy_id%3A11676%29%29+AND+%28reviewed%3Atrue%29'
#
#                     query_uniprot = requests.get(query_string)
#                     query_uniprot = query_uniprot.text.split('\n')
#
#                     query_uniprot_df = pd.DataFrame([x.strip().split('\t') for x in query_uniprot])
#                     cols = query_uniprot_df.iloc[0]
#                     # st.write(cols)
#                     query_uniprot_df = query_uniprot_df[1:len(query_uniprot_df) - 1]
#                     query_uniprot_df.columns = cols
#                     temp = list(query_uniprot_df['Entry'])
#                     # st.write(len(temp))
#                     uprot_list.append(temp)
#
#                 uprot_list = [item for sublist in uprot_list for item in sublist]
#
#                 st.write('A total of', str(len(uprot_list)), 'viral proteins have been identified.', '\n')
#
#                 return (uprot_list)
#

# def createKG(inputs):
#
#     st.write(inputs)
#
#     #if inputs:
#     efo_id = inputs[0]
#     ct_phase = inputs[1]
#     kg_name = inputs[2]
#     st.write(inputs)
#
#     st.write('Now fetching real-time data from databases. Be patient!')
#     st.write('\n')
#
#     temp_id = efo_id.split(' ')
#     # st.write(temp_id)
#     temp_id = [int(x) for x in temp_id]
#     st.write(temp_id)
#     st.write(ct_phase)
#
#     drugs_df = pd.DataFrame()
#     dis2prot_df = pd.DataFrame()
#     dis2snp_df = pd.DataFrame()
#
#     st.write(doid)
#
#     for id in temp_id:
#         chembl_list = GetDiseaseAssociatedDrugs(doid['id'][id], ct_phase)
#         drugs_df = pd.concat([drugs_df, chembl_list])
#
#         prot_list = GetDiseaseAssociatedProteins(doid['id'][id])
#         dis2prot_df = pd.concat([dis2prot_df, prot_list])
#
#         snp_dgnet = GetDiseaseSNPs(doid['id'][id])
#         dis2snp_df = pd.concat([dis2snp_df, snp_dgnet])
#
#     drugs_df = drugs_df.reset_index(drop=True)
#
#     dis2prot_df = dis2prot_df.reset_index(drop=True)
#     dis2snp_df = dis2snp_df.reset_index(drop=True)
#
#     uprot_df = GetDiseaseAssociatedProteinsPlot(dis2prot_df)
#     adv_effect = pd.DataFrame()
#
#     # create empty KG
#     kg = pybel.BELGraph(name=kg_name, version="0.0.1")
#
#     uprot_ext = ExtractFromUniProt(list(set(uprot_df['UniProt'])))
#
#     if vir_prot:
#         vir_uprot_ext = ExtractFromUniProt(vir_prot)
#
#     if not drugs_df.empty:
#         st.write('A total of ' + str(
#             len(list(set(drugs_df['drugId'])))) + ' drugs have been identified. Now fetching relevant data!')
#
#         chembl2mech = RetMech(list(set(drugs_df['drugId'])))
#         chembl2act = RetAct(list(set(drugs_df['drugId'])))
#
#         prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(chembl2mech)
#         prtn_as_chembl = set(prtn_as_chembl)
#         prtn_as_chembl = list(prtn_as_chembl)
#         chembl2uprot = chembl2uniprot(prtn_as_chembl)
#
#         chembl2act = chembl2gene2path(chembl2uprot, chembl2act)
#         chembl2mech = chembl2gene2path(chembl2uprot, chembl2mech)
#
#         kg = chem2moa_rel(chembl2mech, 'HGNC', kg)
#         kg = chem2act_rel(chembl2act, 'HGNC', kg)
#         kg = gene2path_rel(chembl2uprot, 'HGNC', kg)
#
#         adv_effect = GetAdverseEvents(list(set(drugs_df['drugId'])))
#         kg = chembl2adverseEffect_rel(adv_effect, kg)
#
#     kg = uniprot_rel(uprot_ext, 'HGNC', kg)
#
#     if vir_prot:
#         kg = uniprot_rel(vir_uprot_ext, 'VP', kg)
#
#     # snp_dgnet = GetDiseaseSNPs(doid['id'][efo_id])
#
#     # if snp_dgnet != None:
#     # kg = snp2gene_rel(snp_dgnet,kg)
#
#     if not dis2snp_df.empty:
#         kg = snp2gene_rel(dis2snp_df, kg)
#
#     st.write('Your KG is now generated!', '\n')
#
#     #saveFiles(kg_name, uprot_df, adv_effect, kg, drugs_df, dis2snp_df, uprot_ext)
#
#     return (kg)


def GetDiseaseAssociatedProteins(disease_id):

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

    # variables = {"$efo_id":efo_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    # r = requests.post(base_url, json={"query": query_string, "variables": variables})
    r = requests.post(base_url, json={"query": query_string})
    # st.write(r.status_code)

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    temp_list = []
    for item in api_response["data"]["disease"]["associatedTargets"]["rows"]:
        # st.write(item['target'])
        # break
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
    st.write(
        "We have identified "
        + str(len(df))
        + " proteins (Swiss-Prot) associated with the disease. Please note that the proteins identified may not be unique if you combined two or more diseases. Following is a histogram that shows "
        + "distribution of proteins based on scores provided by OpenTargets. The scores are influenced by various factors "
        + "such as genetic associations, expression, mutations, known pathways, targeting drugs and so on."
        + "\n"
    )

    st.write(
        "A total of "
        + str(len(list(set(df["UniProt"]))))
        + " unique proteins have been identified."
    )

    st.write("\n")

    st.write("Displaying top 20 genes")

    fig, ax = plt.subplots()
    ax.hist(df["Score"])
    ax.set_title("Distribution of proteins based on OpenTargets score")
    ax.set_xlabel("Score")
    ax.set_ylabel("No. of proteins")

    fig.tight_layout()
    # st.pyplot(plt.gcf(),use_container_width=True)

    st.header("Summary of proteins related to the disease(s)", divider="gray")

    # title_container = st.beta_container()
    # col = st.columns((1.5, 1.5), gap="medium")
    # #col1, col2 = st.columns(2)
    # #image = Image.open('/home/ddutt/Pictures/Suzieq-logo-2.jpg')
    # #with title_container:
    # with col[0]:
    #     st.write(df.head(20))
    # with col[1]:
    #     st.pyplot(plt.gcf())

    col1, col2 = st.columns(2)
    with col1:
        st.write(df.head(20))
    with col2:
        st.pyplot(fig)

    time.sleep(0.05)

    # with st.form(key='proteinScores'):
    #     score = st.number_input("We recommend taking a threshold above 0.3 to exclude loosely associated proteins. Please enter your desired threshold: ")
    #
    #     submit_button = st.form_submit_button(label='Submit score')

    score = st.number_input(
        "Enter threshold score (recommended > 0.3):",
        min_value=0.0,
        max_value=1.0,
        value=0.3,
        step=0.1,
    )

    if st.button("Apply Score"):
        filtered_df = df[df["Score"] >= score]
        st.success(f"Score applied: {score}. Generating final KG...", icon="✅")
        st.write(f"Total proteins after filtering: {len(filtered_df)}")
        st.session_state.filtered_protein_df = filtered_df
        st.write(filtered_df)
        kg_ph = finalizeKG(filtered_df)
        st.write(kg_ph.summarize())
        return kg_ph

        # st.session_state.kg_finalized = True
        # st.rerun()

    # if submit_button('Apply score'):
    #     st.success('Score saved', icon="✅")
    #     st.write('Score selected: ', score)
    #
    # # score = st.number_input(
    # #     "We recommend taking a threshold above 0.3 to exclude loosely associated proteins. Please enter your desired threshold: ")
    #
    # #st.write('Score selected: ',score)
    #
    #     filter_df = df.loc[df['Score'] >= score, :]
    #
    #     st.write('\n')
    #     st.write('Alright, we are good to go now. Your KG is now being generated! Sit back and relax!!')
    #
    #     st.write('\n', 'Total no. of proteins after filtering with score: ', len(filter_df))
    #
    #     return(filter_df)

    # display(HTML(df.to_html(index=False)))
    # st.write(df)
    # st.write('\n')
    # st.write('ready to return proteinDF')
    # #st.write(df)
    # return (df)


def getDrugCount(disease_id):
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

    # replace $efo_id with value from efo_id
    # query_string = query_string.replace("$efo_id",f'"{efo_id}"')
    # query_string = query_string.replace("$efo_id",f'"{efo_id}"')

    # Set variables object of arguments to be passed to endpoint
    variables = {"my_efo_id": efo_id, "my_size": size}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string, "variables": variables})
    # r = requests.post(base_url, json={"query": query_string})
    # st.write(r.status_code)

    # Transform API response from JSON into Python dictionary and print in console
    api_response = json.loads(r.text)

    df = pd.DataFrame(api_response["data"]["disease"]["knownDrugs"]["rows"])

    # df = df.loc[df['phase'] >= int(CT_phase),:]

    if not df.empty:
        df = df.loc[df["phase"] >= int(CT_phase), :]
        df["id"] = efo_id
        df["disease"] = api_response["data"]["disease"]["name"]
        # st.write('Your dataframe is ready')
        return df

    else:
        st.write("No drugs found in clinical trials")


def KG_namespace_plot(final_kg, kg_name):

    nspace_count = pybel.struct.summary.count_namespaces(final_kg)
    nspace_count = dict(nspace_count)

    nspace_data = {
        "Namespace": list(nspace_count.keys()),
        "Number": list(nspace_count.values()),
    }
    nspace = pd.DataFrame(nspace_data)
    plt.figure()

    a = sns.barplot(x="Number", y="Namespace", data=nspace_data)
    a.set(xlabel="Number", ylabel="Namespace", title="KG Namespace in numbers")

    plt.tight_layout()
    plt.savefig(kg_name + "_namespace.png", dpi=600)
    # plt.show()


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
            # st.write(r.status_code)

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
        graph.add_association(
            Abundance(namespace="ChEMBL", name=str(chembl_adveff_df["chembl_id"][i])),
            Pathology(
                namespace="SideEffect", name=str(chembl_adveff_df["name"][i])
            ),  # TODO: Fix namespace
            citation="OpenTargets Platform",
            evidence="DrugReactions",
        )

    return graph


def GetViralProteins(query_disease):
    # file downloaded from https://www.genome.jp/ftp/db/virushostdb Dated: 12/09/2023
    virus = pd.read_csv(
        "https://raw.githubusercontent.com/Fraunhofer-ITMP/kgg/main/data/virushostdb.csv"
    )

    cols = ["virus tax id", "virus name", "DISEASE", "host tax id"]
    virus = virus[cols]

    # st.write(virus)

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

        # st.write(disease)
        st.write("\n")
        st.write(
            "The workflow has identified your query as a viral disease. Its proteins (SWISS-Prot) will be now represented in the KG.",
            "\n",
        )

        st.write("Enter Viral Prot Function")

        time.sleep(0.1)

        with st.form(key="identify viral protein"):

            virus_name = st.text_input(
                "Do you want to look further for a specific virus? Please type its name or skip it by typing 'no': "
            )

            submit_button = st.form_submit_button(label="Search")

        if submit_button:

            # virus_name = virus_name

            # virus_name = input(
            #     'Do you want to look further for a specific virus? Please type its name or skip it by typing \'no\': ')

            # subset df with virus name
            if virus_name.lower() != "no":

                virus_subset_2 = virus[
                    virus["virus name"].str.contains(virus_name, na=False, case=False)
                ]
                # st.write(virus_subset_2)
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
            st.write(virus_subset_merge)

            # display(HTML(virus_subset_merge.to_html(index=False)))

            time.sleep(0.1)

            temp_id = st.text_input(
                "Enter the index value(s). If multiple, use space, for example -> 0 1 3: "
            )

            # temp_id = input('Enter the index value(s). If multiple, use space, for example -> 0 1 3: ')

            st.write("\n")

            temp_id = temp_id.split(" ")
            temp_id = [int(x) for x in temp_id]
            # st.write(virus_subset_merge.loc[0]['virus tax id'])

            uprot_list = []

            for item in temp_id:
                tax_id = virus_subset_merge.loc[item]["virus tax id"]
                # st.write(tax_id)

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
                # st.write(cols)
                query_uniprot_df = query_uniprot_df[1 : len(query_uniprot_df) - 1]
                query_uniprot_df.columns = cols
                temp = list(query_uniprot_df["Entry"])
                # st.write(len(temp))
                uprot_list.append(temp)

            uprot_list = [item for sublist in uprot_list for item in sublist]

            st.write(
                "A total of",
                str(len(uprot_list)),
                "viral proteins have been identified.",
                "\n",
            )

            return uprot_list


def saveFiles(
    kgName, disease2protein, drugAdvEffect, final_kg, drug_df, snp_df, uprot_dict
):
    st.write("Now let's save all the files that were created in the process.", "\n")

    st.write(
        "Please enter the location (e.g. 'C:\\Users\\rkarki\\Documents\\kg\\' ) where KG files should be stored. A folder will be created automatically.",
        "\n",
    )

    time.sleep(0.2)

    path = st.text_input("Your desired location")

    path = path + kgName

    original_path = sys.path[0]

    os.makedirs(path, exist_ok=True)

    os.chdir(path)

    disease2protein.to_csv("diseaseAssociatedProteins.csv", sep=",")

    drugAdvEffect.to_csv("adverseEffects.csv", sep=",")

    drug_df.to_csv("diseaseAssociatedDrugs.csv", sep=",")

    snp_df.to_csv("diseaseAssociatedSNPs.csv", sep=",")

    # to cytoscape compatible graphml
    pybel.to_graphml(final_kg, kgName + ".graphml")

    # to regular BEL format
    pybel.dump(final_kg, kgName + ".bel")

    # to neo4j
    pybel.to_csv(final_kg, kgName + ".csv")

    filename = kgName + ".pkl"
    outfile = open(filename, "wb")
    pickle.dump(final_kg, outfile)
    outfile.close()

    filename = kgName + "_prot_dict" + ".pkl"
    outfile = open(filename, "wb")
    pickle.dump(uprot_dict, outfile)
    outfile.close()

    # plot for namespace dist
    KG_namespace_plot(final_kg, kgName)

    os.chdir(original_path)


# @animation.wait(wheel,text = '    Fetching SNPs ')
def GetDiseaseSNPs(disease_id):
    try:
        snps = get_variants_by_efo_id(disease_id)
        snps_df = snps.genomic_contexts
        snps_df["disease_id"] = disease_id
        # snps_functional_class_df = snps.variants
        # snps_functional_class_df['disease_id'] = disease_id

        return snps_df

    except:
        st.write(
            "No SNPs found. This could be either because 1. The identifier of your disease of interest is not compatible with Experimental Factor Ontology (EFO) or 2. No SNPs have been reported for the disease"
        )


def snp2gene_rel(snp_df, graph):
    # st.write(snp_df.head(2))

    kg_prots = getProtfromKG(graph)

    unique_prots_df = pd.DataFrame(kg_prots, columns=["Proteins"])

    snp_df = unique_prots_df.merge(
        snp_df, how="inner", left_on="Proteins", right_on="gene.geneName"
    )

    # take SNPs which are within the gene sequence
    snp_df = snp_df.loc[snp_df["distance"] == 0]
    snp_df = snp_df.reset_index(drop=True)

    st.write(
        "A total of "
        + str(len(snp_df))
        + " SNPs have been identified from GWAS Central. Now adding relevant data"
    )
    st.write("\n")

    # graph = pybel.BELGraph(name='test', version="0.0.1")

    for i in stqdm(range(len(snp_df)), desc="Adding disease associated SNPs"):
        graph.add_association(
            Gene(namespace="dbSNP", name=snp_df["rsId"][i]),
            Protein(namespace="HGNC", name=snp_df["Proteins"][i]),
            citation="GWAS Central",
            evidence="SNPs for queried disease",
        )

    return graph


def getViralProteins(query_disease):
    # File loading and initial processing (unchanged)
    virus = pd.read_csv(
        "https://raw.githubusercontent.com/Fraunhofer-ITMP/kgg/main/data/virushostdb.csv"
    )
    cols = ["virus tax id", "virus name", "DISEASE", "host tax id"]
    virus = virus[cols]
    virus = virus.loc[virus["host tax id"] == 9606.0, :]
    virus = virus.reset_index(drop=True)
    virus["host tax id"] = pd.to_numeric(virus["host tax id"], downcast="integer")

    # Subset df with disease keyword
    virus_subset_1 = virus[
        virus["DISEASE"].str.contains(query_disease, na=False, case=False)
    ]

    if not virus_subset_1.empty:
        st.write("\n")
        st.write(
            "The workflow has identified your query as a viral disease. Its proteins (SWISS-Prot) will be now represented in the KG.",
            "\n",
        )
        st.write("Enter Viral Prot Function")

        with st.form(key="identify_viral_proteins"):
            virus_name = st.text_input(
                "Do you want to look further for a specific virus? Please type its name or skip it by typing 'no': ",
                value=st.session_state.virus_name,
            )
            submit_button = st.form_submit_button(
                label="Search", on_click=virus_form_callback()
            )

        if st.session_state.virus_form_submitted:
            st.write(st.session_state.virus_name, " submit button")

            if st.session_state.virus_name.lower() != "no":
                st.session_state.virus_subset_2 = virus[
                    virus["virus name"].str.contains(
                        st.session_state.virus_name, na=False, case=False
                    )
                ]
                if not st.session_state.virus_subset_2.empty:
                    st.write("Virus subset based on your input:")
                    st.write(st.session_state.virus_subset_2)
                else:
                    st.write("No matching viruses found for the given name.")
            else:
                st.write("Showing all results for the disease:")
                st.write(virus_subset_1)

            # Reset the form submitted flag
            st.session_state.virus_form_submitted = False

    return virus_subset_1


def createInitialKG(inputs):
    st.write(inputs)

    efo_id = inputs["efo_idx"]
    ct_phase = inputs["ct_phase"]
    kg_name = inputs["namekg"]
    doid = inputs["df"]

    st.write("Now fetching real-time data from databases. Be patient!")

    temp_id = efo_id.split(" ")
    temp_id = [int(x) for x in temp_id]

    drugs_df = pd.DataFrame()
    dis2prot_df = pd.DataFrame()
    dis2snp_df = pd.DataFrame()

    for id in temp_id:
        st.write("enter for loop")

        chembl_list = GetDiseaseAssociatedDrugs(doid["id"][id], ct_phase)
        drugs_df = pd.concat([drugs_df, chembl_list])

        prot_list = GetDiseaseAssociatedProteins(doid["id"][id])
        dis2prot_df = pd.concat([dis2prot_df, prot_list])

        snp_dgnet = GetDiseaseSNPs(doid["id"][id])
        dis2snp_df = pd.concat([dis2snp_df, snp_dgnet])

    st.write("exit for loop")

    drugs_df = drugs_df.reset_index(drop=True)
    dis2prot_df = dis2prot_df.reset_index(drop=True)
    dis2snp_df = dis2snp_df.reset_index(drop=True)

    st.session_state.drugs_df = drugs_df
    st.session_state.dis2prot_df = dis2prot_df
    st.session_state.dis2snp_df = dis2snp_df
    st.session_state.kg_name = kg_name
    st.session_state.initial_kg_created = True


def finalizeKG(uprot_df):
    # This function will contain the remaining code from your original createKG function
    # that needs to run after protein filtering

    kg_name = st.session_state.kg_name
    vir_prot = None
    drugs_df = st.session_state.drugs_df
    st.write(drugs_df)

    st.write(len(drugs_df))

    dis2snp_df = st.session_state.dis2snp_df

    st.write("Finalizing KG...")
    adv_effect = pd.DataFrame()

    # create empty KG
    kg = pybel.BELGraph(name=kg_name, version="0.0.1")

    uprot_ext = ExtractFromUniProt(list(set(uprot_df["UniProt"])))
    # st.write(uprot_ext)

    if vir_prot:
        vir_uprot_ext = ExtractFromUniProt(vir_prot)

    if not drugs_df.empty:
        st.write(
            "A total of "
            + str(len(list(set(drugs_df["drugId"]))))
            + " drugs have been identified. Now fetching relevant data!"
        )

        chembl2mech = RetMech(list(set(drugs_df["drugId"])))
        chembl2act = RetAct(list(set(drugs_df["drugId"])))

        prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(
            chembl2mech
        )
        prtn_as_chembl = set(prtn_as_chembl)
        prtn_as_chembl = list(prtn_as_chembl)
        chembl2uprot = chembl2uniprot(prtn_as_chembl)

        chembl2act = chembl2gene2path(chembl2uprot, chembl2act)
        chembl2mech = chembl2gene2path(chembl2uprot, chembl2mech)

        kg = chem2moa_rel(chembl2mech, "HGNC", kg)
        kg = chem2act_rel(chembl2act, "HGNC", kg)
        kg = gene2path_rel(chembl2uprot, "HGNC", kg)

        adv_effect = GetAdverseEvents(list(set(drugs_df["drugId"])))
        kg = chembl2adverseEffect_rel(adv_effect, kg)

    kg = uniprot_rel(uprot_ext, "HGNC", kg)

    if vir_prot:
        kg = uniprot_rel(vir_uprot_ext, "VP", kg)

    # snp_dgnet = GetDiseaseSNPs(doid['id'][efo_id])

    # if snp_dgnet != None:
    # kg = snp2gene_rel(snp_dgnet,kg)

    if not dis2snp_df.empty:
        kg = snp2gene_rel(dis2snp_df, kg)

    st.write("Your KG is now generated!", "\n")

    # saveFiles(kg_name, uprot_df, adv_effect, kg, drugs_df, dis2snp_df, uprot_ext)

    # return (kg)
    # Your remaining code here
    # You can access st.session_state.filtered_protein_df, st.session_state.drugs_df, etc.

    st.write("KG creation completed!")
    st.session_state.kg_finalized = True
    st.write("KG summary: ", kg.summarize())
    return kg


def createKG(inputs):

    st.write(inputs)

    # if inputs:
    efo_id = inputs["efo_idx"]
    ct_phase = inputs["ct_phase"]
    kg_name = inputs["namekg"]
    doid = inputs["df"]

    st.write("Now fetching real-time data from databases. Be patient!")
    st.write("\n")

    temp_id = efo_id.split(" ")
    # st.write(temp_id)
    temp_id = [int(x) for x in temp_id]

    st.write(temp_id)
    st.write(ct_phase)

    drugs_df = pd.DataFrame()
    dis2prot_df = pd.DataFrame()
    dis2snp_df = pd.DataFrame()

    st.write(doid)

    for id in temp_id:
        st.write("enter for loop")

        chembl_list = GetDiseaseAssociatedDrugs(doid["id"][id], ct_phase)
        drugs_df = pd.concat([drugs_df, chembl_list])

        prot_list = GetDiseaseAssociatedProteins(doid["id"][id])
        dis2prot_df = pd.concat([dis2prot_df, prot_list])

        snp_dgnet = GetDiseaseSNPs(doid["id"][id])
        dis2snp_df = pd.concat([dis2snp_df, snp_dgnet])

    st.write("exit for loop")

    drugs_df = drugs_df.reset_index(drop=True)

    dis2prot_df = dis2prot_df.reset_index(drop=True)
    dis2snp_df = dis2snp_df.reset_index(drop=True)

    st.session_state.drugs_df = drugs_df
    st.session_state.dis2prot_df = dis2prot_df
    st.session_state.dis2snp_df = dis2snp_df
    st.session_state.kg_created = True

    uprot_df = GetDiseaseAssociatedProteinsPlot(dis2prot_df)
    adv_effect = pd.DataFrame()

    # create empty KG
    kg = pybel.BELGraph(name=kg_name, version="0.0.1")

    uprot_ext = ExtractFromUniProt(list(set(uprot_df["UniProt"])))

    if vir_prot:
        vir_uprot_ext = ExtractFromUniProt(vir_prot)

    if not drugs_df.empty:
        st.write(
            "A total of "
            + str(len(list(set(drugs_df["drugId"]))))
            + " drugs have been identified. Now fetching relevant data!"
        )

        chembl2mech = RetMech(list(set(drugs_df["drugId"])))
        chembl2act = RetAct(list(set(drugs_df["drugId"])))

        prtn_as_chembl = Ret_chembl_protein(chembl2act) + Ret_chembl_protein(
            chembl2mech
        )
        prtn_as_chembl = set(prtn_as_chembl)
        prtn_as_chembl = list(prtn_as_chembl)
        chembl2uprot = chembl2uniprot(prtn_as_chembl)

        chembl2act = chembl2gene2path(chembl2uprot, chembl2act)
        chembl2mech = chembl2gene2path(chembl2uprot, chembl2mech)

        kg = chem2moa_rel(chembl2mech, "HGNC", kg)
        kg = chem2act_rel(chembl2act, "HGNC", kg)
        kg = gene2path_rel(chembl2uprot, "HGNC", kg)

        adv_effect = GetAdverseEvents(list(set(drugs_df["drugId"])))
        kg = chembl2adverseEffect_rel(adv_effect, kg)

    kg = uniprot_rel(uprot_ext, "HGNC", kg)

    if vir_prot:
        kg = uniprot_rel(vir_uprot_ext, "VP", kg)

    # snp_dgnet = GetDiseaseSNPs(doid['id'][efo_id])

    # if snp_dgnet != None:
    # kg = snp2gene_rel(snp_dgnet,kg)

    if not dis2snp_df.empty:
        kg = snp2gene_rel(dis2snp_df, kg)

    st.write("Your KG is now generated!", "\n")

    # saveFiles(kg_name, uprot_df, adv_effect, kg, drugs_df, dis2snp_df, uprot_ext)

    return kg


# Initialize session state variables
if "query_generated" not in st.session_state:
    st.session_state.query_generated = False
    st.session_state.query = None
    # st.session_state.disease = None
    st.session_state.search_df = None

# Add a button to clear the previous search and start a new one
if st.session_state.query_generated:
    if st.button("Search for another disease"):
        st.session_state.query_generated = False
        st.session_state.query = None
        # st.session_state.disease = None
        st.session_state.search_df = None
#
# Only call Generate_KG() if kg_generated is False
if not st.session_state.query_generated:
    query = getQuery()
    if query:
        disease_df(query)
#
# Check if kg_generated is True before using the results
if st.session_state.query_generated:
    # Process the results

    st.write("let/'s do something more")
    st.write(st.session_state.query)
    get_inputs()
    # st.write('Printing', get_inputs())

    # st.write(st.session_state.query)
    # st.write(st.session_state.search_df)

if "query_phaseCompleted" not in st.session_state:
    st.session_state.query_phaseCompleted = False
    st.session_state.efo_id = None
    st.session_state.ctphase = None
    st.session_state.kg_name = None
#
if st.session_state.query_phaseCompleted:
    if st.button("Change inputs"):

        del st.session_state.efo_id
        del st.session_state.ctphase
        del st.session_state.kg_name
        del st.session_state.query_phaseCompleted

        st.session_state.query_phaseCompleted = False
        st.session_state.efo_id = None
        st.session_state.ctphase = None
        st.session_state.kg_name = None

# if st.session_state.query_phaseCompleted:
#
#     #st.write(st.session_state)
#
#     if st.button('Create KG'):
#
#         inputs = {'query':st.session_state.query, 'df':st.session_state.search_df,
#         'efo_idx':st.session_state.efo_id,'ct_phase':st.session_state.ctphase,'namekg':st.session_state.kg_name}
#
#         st.write(len(inputs))
#
#         createKG(inputs)


# # Main app logic
# if 'kg_created' not in st.session_state:
#     st.session_state.kg_created = False
#
# if 'kg_finalized' not in st.session_state:
#     st.session_state.kg_finalized = False
#
# if st.session_state.query_phaseCompleted and not st.session_state.kg_created:
#     if st.button('Create KG'):
#         inputs = {
#             'query': st.session_state.query,
#             'df': st.session_state.search_df,
#             'efo_idx': st.session_state.efo_id,
#             'ct_phase': st.session_state.ctphase,
#             'namekg': st.session_state.kg_name
#         }
#         createKG(inputs)
#
# if st.session_state.kg_created and not st.session_state.kg_finalized:
#     GetDiseaseAssociatedProteinsPlot(st.session_state.dis2prot_df)
#
# if st.session_state.kg_finalized:
#     st.write("Final KG created with filtered proteins:")
#     st.write(st.session_state.filtered_protein_df)
#
#     # Add option to start over
#     if st.button('Start New KG'):
#         for key in ['query_generated', 'query', 'search_df', 'query_phaseCompleted', 'efo_id', 'ctphase', 'kg_name',
#                     'kg_created', 'kg_finalized', 'drugs_df', 'dis2prot_df', 'dis2snp_df', 'filtered_protein_df']:
#             if key in st.session_state:
#                 del st.session_state[key]
#         st.rerun()


# Main app logic
if "initial_kg_created" not in st.session_state:
    st.session_state.initial_kg_created = False

if "kg_finalized" not in st.session_state:
    st.session_state.kg_finalized = False

if st.session_state.query_phaseCompleted and not st.session_state.initial_kg_created:
    if st.button("Create Initial KG"):
        inputs = {
            "query": st.session_state.query,
            "df": st.session_state.search_df,
            "efo_idx": st.session_state.efo_id,
            "ct_phase": st.session_state.ctphase,
            "namekg": st.session_state.kg_name,
        }
        createInitialKG(inputs)
        st.rerun()

if st.session_state.initial_kg_created and not st.session_state.kg_finalized:
    GetDiseaseAssociatedProteinsPlot(st.session_state.dis2prot_df)

if st.session_state.kg_finalized:
    st.write("KG creation process completed!")
    # Display final KG information or provide download options

    # Add option to start over
    if st.button("Start New KG"):
        for key in [
            "query_generated",
            "query",
            "search_df",
            "query_phaseCompleted",
            "efo_id",
            "ctphase",
            "kg_name",
            "initial_kg_created",
            "kg_finalized",
            "drugs_df",
            "dis2prot_df",
            "dis2snp_df",
        ]:
            if key in st.session_state:
                del st.session_state[key]
        st.rerun()
