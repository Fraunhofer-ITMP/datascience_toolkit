#!/bin/bash

uvicorn kgg_api:app --host 0.0.0.0 --port 8080 --log-level info --reload