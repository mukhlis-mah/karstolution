try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description':"""Karstolution module: a speleothem d18O Proxy System Module (PSM) 
    that incorporates Karst and in-cave isotope fractionation processes (coupling of KarstFor and ISOLUTION). 
    Windows executable version with a GUI also avaiable.""",
    'author':"Mukhlis Mah",
    'url':'https://github.com/mukhlis-mah/karstolution',
    'download_url':'https://github.com/mukhlis-mah/karstolution',
    'author_email':'mukhlis.mah@gmail.com',
    'version':'0.1',
    'install_requires':['nose','numpy','scipy'],
    'packages':['Karstolution'],
    'scripts':['cmodel_frac.py','constants.py','evaporation.py','isotope_calcite.py',
    'karst_process.py','karstolution1.1','O18EVA.py','O18EVA_MEAN.py'],
    'name':'Karstolution'
}

setup(**config)