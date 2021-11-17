# Intrinsically disordered protein simulation database

This simulation database library is used to manage and analyse the large SolSpace intriniscally disorderded protein (IDP) database.

In this library, each protein will be simulated under 3-9 different solution condition. Each solution condition will
have 3-5 repeats. I created a class to automatically detect the simulation configuration and directory of the dataset.

This library can also generate a dataframe/csv file to summarize the structural property of IDPs including my novel
interaction map method.

# What is intriniscally disordered protein?

(Brief introduction of IDP. Copy from abstract)

I simulated over 200 IDP sequences based on our SolSpace all-atom simulation. 

# My goal:

1. Discover how surrounding environment will alter IDP structure.
2. Predict IDP solution sensitivity and structual change based on its sequence
3. Design single point mutation based my prediction to expand or compact the IDP conformational ensemble.
