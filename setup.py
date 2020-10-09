from setuptools import setup, find_packages

setup(
    name='IS2_velocity',
    version='0.0.1',
    license='BSD',
    description='surface ice velocity mapping with ICESat-2',
    url='https://github.com/ICESAT-2HackWeek/IS2_velocity',
    packages=find_packages(exclude=['data','test','notebooks']),
    install_requires=[
        'numpy',
        'scipy',
        'h5py',
        'astropy',
        'icepyx',
        'pointCollection @ git+https://github.com/smithB/pointCollection.git#egg=pointCollection-1.0.0.0'
    ]
)
