# -*- coding: utf-8 -*-

import os
import json
import pandas as pd
import requests
from pybel import BELGraph
from pybel.dsl import Protein, Abundance, Pathology, Gene


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

    return df
