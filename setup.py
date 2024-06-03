from setuptools import setup, find_packages

setup(
    name='magic_irma',
    version='0.1',
    install_requires=[
        'biopython',
        'pandas',
    ],
    entry_points={
        'console_scripts': [
            'irmagic=magic_irma.irmagic:main',
        ],
    },
    author='Iván Bloise Sánchez',
    author_email='ibloise62@gmail.com',
    description='Tools for life-quality working with CDC-IRMA',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/ibloise/magic_irma',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Development Status :: 1 - Planning'
    ],
    python_requires='>=3.6',
)
