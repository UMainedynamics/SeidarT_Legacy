#!/usr/bin/python3

from setuptools import setup

form os import path 
this_directory = path.abspath(path.dirname(__file__) )
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
	long_description = f.read()

setup(
	name = "SeidarT",
	version = "0.0.1",
	author = "Steven Bernsen",
	author_email = "stevenbernsen@gmail.com",
	description = "A combined seismc and radar modeling toolbox",
	long_decscription=long_description,
	long_description_content_type='text/markdown'
)