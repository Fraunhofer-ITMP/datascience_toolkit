# DataScience Toolkit

This repository contains the overall suite of resources and tools built by the Fraunhofer ITMP data science team.

# Local testing

To run the app locally for debuging and testing purposes, do as follows:
```bash
git clone https://github.com/Fraunhofer-ITMP/datascience_toolkit.git
cd datasicence_toolkit
conda create --name=streamlit-app python=3.11
conda activate streamlit-app
pip install -r requirements.txt
streamlit run Main.py
```

This will then open up the streamlit app locally.

# Updating information

To add new members to the group, check out the [data](https://github.com/Fraunhofer-ITMP/datascience_toolkit/tree/main/data) folder. 

* New members: Add the details of the new member in the [CSV](https://github.com/Fraunhofer-ITMP/datascience_toolkit/blob/main/data/members.csv) file.

* New publication: Add a new publication in Chicago style from [Google scholar](https://scholar.google.com/) in the [CSV](https://github.com/Fraunhofer-ITMP/datascience_toolkit/blob/main/data/publications.tsv)

# Updating website

To update the website, you need access to the [SERVE platform](https://serve.scilifelab.se/accounts/login/). If you do not have an account, please create one and reach out to Yojana to assist you in the account creation.

Once you have an account, ask Yojana or Andrea to add you to the project.