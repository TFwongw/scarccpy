from setuptools import setup, find_packages

setup(
    name='scarcc',
    version='0.0.10',
    url='https://github.com/TFwongw/scarccpy',
    author='Thomas Wong',
    author_email='wong0755@umn.edu',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(where="src"),
    package_dir = {"": "src"},
    install_requires=[
        'cobra',
        'pandas<2',
    ],
)