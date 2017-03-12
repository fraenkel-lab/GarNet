from distutils.core import setup

setup(
    name='GarNet',
    version='0.2.1',
    packages=['GarNet'],
    package_dir={'GarNet': 'src'},
    scripts=['src/garnet.py'],
    url='https://github.com/fraenkel-lab/GarNet',
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'],
    license='GNU General Public License',
    author='zfrenchee',
    author_email='alex@lenail.org',
    description='',
    install_requires=[
        "numpy",
        "pandas",
        "statsmodels",
        "matplotlib",
        "intervaltree",
        "jinja2"
    ],
)
