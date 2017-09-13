from setuptools import setup

setup(
    name='GarNet',
    packages=['GarNet'],
    package_data={'GarNet': ['summary.jinja']},
    version='0.5.0',
    url='https://github.com/fraenkel-lab/GarNet',
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'],
    license='MIT',
    author='zfrenchee, iamjli',
    author_email='alex@lenail.org, iamjli@mit.edu',
    description='',
    install_requires=[
        'numpy',
        'pandas',
        'statsmodels',
        'matplotlib',
        'jinja2',
        'pybedtools'
    ],
)

