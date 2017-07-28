from setuptools import setup

setup(
    name='GarNet',
    packages=['GarNet'],
    package_data={'GarNet': ['summary.jinja']},
    version='0.4.5',
    url='https://github.com/fraenkel-lab/GarNet',
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'],
    license='MIT',
    author='zfrenchee',
    author_email='alex@lenail.org',
    description='',
    install_requires=[
        'numpy',
        'pandas',
        'statsmodels',
        'matplotlib',
        'intervaltree',
        'jinja2'
    ],
)

