# QMRA of Photorepair After UV Disinfection

# Attribution:
This QMRA code for photorepair was developed by Daniel Ma with the advising of Dr. Mark H. Weir and Dr. Natalie M. Hull at The Ohio State University. This project was developed as a course project for PUBHES 7375 QMRA Modeling instructed by Dr. Mark H. Weir in Spring 2020 and incorporated into Daniel Ma's MS Thesis in Environmental Engineering at The Ohio State University. 

# Project Dates: 
Feb 2020 - Mar 2023

# Title of Manuscript: 
Fluence-based QMRA model for bacterial photorepair and regrowth in drinking water after decentralized UV disinfection; DOI: https://doi.org/10.1016/j.watres.2023.119612

# What is QMRA?
Check out the quantitative microbial risk assessment wiki (http://qmrawiki.org/) for some more information.

# Project Description:
What this code does: This code calculates the risk of infection following photorepair exposure of drinking water treated by low pressure UV disinfection. The photorepair light exposure is modelled in plastic containers that are carried from the point of collection (e.g. a kiosk, communal water access point, store) to the point of consumption. The collection is modelled as different times of the day and the exposure duration is also modelled as various lengths of time. Different sunlight intensities are accounted for to provide different sunlight scenarios (e.g. times of the year). Different levels of UV disinfection are also modelled and simulated. Photorepair is modelled using a fluence-based model determined from Bohrerova and Linden (2007).

A unique aspect of this code was the need to calculate the photorepair fluence. We based this calculation on spectral (i.e. wavelength dependent) characteristics of light transmission and photorepair enzyme activity. Two sources of solar data were obtained: total irradiance and spectral irradiance, which are described in more detail in the data section below.

Motivation: The motivation of the creating this code was to conduct a formal quantitative microbial risk assessment (QMRA) for photorepair after UV disinfection for packaged water (e.g. water bottling in communal water systems). Photorepair has been studied extensively in the lab. No risk assessments have been performed to quantify risk due to photorepair. A QMRA was performed as a modeling approach to estimate risk in a population.

# Data
Low pressure fluence response data: Log inactivation vs. UV fluence (or UV dose) for various strains of E. coli.
Solar irradiance: Data in the form of W/m2 at various time intervals throughout the day were obtained from website related to solar photovoltaics (https://www.pveducation.org/pvcdrom/properties-of-sunlight/calculation-of-solar-insolation).
Spectral solar irradiance: Data in the form of W/m2/nm - ASTM G-173-03 (International standard ISO 9845-1, 1992).
Material spectra: Transmission values at various wavelengths obtained from the literature and unpublished data from Dr. Natalie Hull.
Photorepair fluence response data: Log reactivation vs. photorepair fluence for E. coli.
Dose response data: Data for E. coli EHEC, dose response parameter (http://qmrawiki.canr.msu.edu/index.php/Escherichia_coli_enterohemorrhagic_(EHEC):_Dose_Response_Models#cite_note-Auld_et_al._2004-2).
