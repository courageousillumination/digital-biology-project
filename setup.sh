#!/bin/bash
virtualenv env
source env
pip install -r requirements.txt
cd data && tar -xvzf phospho.tar.gz
