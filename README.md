# ITSxpress-qiime2
A plugin for qiime2 that runs ITSxpress




# Installation
## Requirments

In order to consider install ITSxpress-qiime2, it is required to install qiime2.

Here is the installation page to [install qiime](https://docs.qiime2.org/2018.6/install/)

Once qiime2 has been install, you can now install the dependices for ITSxpress.

1. Add the bioconda channel to the qiime2.
  
		conda config --add channels bioconda
			 
2. Install hmmer.
	
		conda install hmmer
		
3. Install bbmap.

		conda install bbmap
	
4. Install ITSxpress

		conda instal bbmap
## ITSxpress-qiime2

To install the ITSxpress plugin for qiime, a few steps are needed.

1. Clone the project onto your computer.

		git clone https://github.com/kweber1/ITSxpress-qiime2.git
		
2. Install the project using pip.

		pip install .

		pip install -e .
		
3. Open your qiime2 environment.
	
		qiime dev refresh-cache
		
4. Check to see if the ITSxpress plugin is installed.

		qiime itsxpress
![Image of console](https://i.gyazo.com/bc013672a324123209b284f889eaa277.png)






