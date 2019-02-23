# Reliable Crop Identification with Satellite Imagery in the Context of Common Agriculture Policy Subsidy Control

R code used for Landsat image proccessing, KNN and SVM classifier training and validation, and calibration method development.

Code documented in part English and part Portuguese. Translations might be available at request.

### About the paper

**Authors:** Jonas Schmedtmann and Manuel Campagnolo, School of Agriculture, University of Lisbon, Portugal.

Full-text available at [https://www.mdpi.com/2072-4292/7/7/9325](https://www.mdpi.com/2072-4292/7/7/9325).

Full Master's Thesis available at [ResearchGate](https://www.researchgate.net/publication/269333638_Automatizing_photo_interpretation_of_satellite_imagery_in_the_context_of_the_Common_Agriculture_Policy_subsidy_control).

### Abstract

Agricultural subsidies in the context of the Common Agricultural Policy (CAP) represent over 40% of the EU’s yearly budget. To ensure that funds are properly spent, farmers are controlled by National Control and Paying Agencies (NCPA) using tools, such as computer-assisted photo interpretation (CAPI), which aims at identifying crops via remotely-sensed imagery. CAPI is time consuming and requires a large team of skilled photo interpreters. The objective of this study was to develop a reliable control system to partially replace CAPI for crop identification, with the overreaching goal of reducing control costs and completion time. Validated control data provided by the Portuguese Control and Paying Agency and an atmospherically-corrected Landsat ETM+ time series were used to perform parcel-based crop classification, leading to an accuracy of only 68% due to high similarity between crops’ spectral signatures. To address this problem, we propose an automatic control system (ACS) that couples crop classification to a reliability requirement. This allows the decision-maker to set a reliability level, which restricts automatic crop identification to parcels that are classified with high certainty. While higher reliability levels reduce the risk of misclassifications, lower levels increase the proportion of automatic control decisions (ACP). With a reliability level of 80%, more than half of the parcels in our study area are automatically identified with an overall accuracy of 84%. In particular, this allows automatically controlling over 85% of all parcels classified as maize, rice, wheat or vineyard.
