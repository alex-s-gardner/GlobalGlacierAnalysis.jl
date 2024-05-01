using FTPClient

# Cite As: "EOLIS elevation data generated using swath processing of CryoSat data (Gourmelen, N., Escorihuela, M., Shepherd, A., Foresta, L., Muir, A., Garcia-Mondejar, A., Roca, M., Baker, # # S., & Drinkwater, M. R. (2018)) and provided by the ESA CryoTEMPO project (https://cryotempo-eolis.org/)."

# chmod 600 ~/.netrc
#;aria2c ftp://science-pds.cryosat.esa.int//TEMPO_SWATH_POINT -c -d /mnt/bylot-r3/data/


#aria2c ftp://science-pds.cryosat.esa.int/TEMPO_SWATH_POINT/2010/07/GREENLAND/CS_OFFL_THEM_POINT_GREENLAND_2010_07_-100000___-800000___V002.nc -c -d /mnt/bylot-r3/data/

ftp = FTP(hostname = "science-pds.cryosat.esa.int/TEMPO_SWATH_POINT/2010/07/", username = "anonymous", password = "nothing")
r = readdir(ftp)

cd(ftp, "GREENLAND")
