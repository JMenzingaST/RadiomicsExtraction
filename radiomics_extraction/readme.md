Within this folder various scripts work together to load dicom files and series, extract information out of them and perform a number of manipulations so that in the end, you end up with Radiomics features.

These features can then in turn be used in a ML project, like for instance, improving Normal Tissue Complication Probabilty models.

I left a few variations of main scripts in here. During data collection it was practical to create seperate scripts to collect data from a certain subsection. I hope that leaving them as examples also serves to show any future users the possibillities of te code provided.

The below image shows how the different components work together:
![Flowchart depicting how the different classes work together to extract radiomics features from a CT scan]([radiomics_extraction\Script_flowchart.png](https://github.com/JMenzingaST/RadiomicsExtraction/blob/main/radiomics_extraction/Script_flowchart.png)?raw=true "Script Flowchart")

Creating this pipeline took a lot of time, trial and error, frustration and perseverance to get working. I hope it will serve future researchers well.

