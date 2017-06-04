Serial version of the program, where the processing of the data is NOT made during the simulation. 
Vacancies are in seperate array. This implementation does not use the SPRNG library.

--------------------------------------------------------------------------------------------------
This version will make 3 subfolders and one file

siminfo.txt - information about the simulation 
data - where the CA is dumped
img - where .ppm images are generated
postdata - where data about the clusters will be generated (e.g. step shortclusdata.txt)
--------------------------------------------------------------------------------------------------
To compile the programs use the bash script compile.sh

./compile.sh

To run the programs first

./npv.x

then

./post.sh
