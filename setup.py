from setuptools import setup, find_packages

setup(
    name='scarcc',
    version='0.1.0',
    url='https://github.com/TFwongw/scarccpy',
    author='Thomas Wong',
    author_email='wong0755@umn.edu',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(where="src"),
    package_dir = {"": "src"},
    python_requires='>=3.10',
    install_requires=[
        'pandas<2.5',
        'matplotlib',
        'seaborn>0.11',
        'cobra',
        'cometspy'
    ],
)