"""
Module containing information to do with chemical speices.
Allows you to obtain latex stings and RMM or smiles from a GeosChem species.
"""

class species:
    """
    Class containing various informaiton about a chemical species, including its RMM
    and latex string.
    """

    def __repr__(self):
         return "This is a class to hold chemical species information"
    def __str__(self):
        try: 
            print "name = " + self.name
        except:
            pass
        try: 
            print "GeosChem group = " + self.group
        except:
            pass
        try:
             print "formula = " + self.formula
        except:
            pass
        try:
             print "InchI = " + self.InChI
        except:
            pass
        try:
             print "smiles = " + self.smiles
        except:
            pass
        try:
             print "RMM = " + str(self.RMM)
        except:
            pass
        try:
             print "Latex = " + self.Latex
        except:
            pass
        return "Class for holding species data"

    def __init__(self, name):
        """
        Function to initilise the class.
        """
        import os
        import csv
        self.name = name
        species_filename = os.path.dirname(__file__) + "/data/species.csv"

        try:
            species_file     = open(species_filename, 'rb')
        except IOError:
            print "Error: Species.csv does not appear to exist."
        species_csv      = csv.reader(species_file)

        if (name in ['OH', 'HO2']):
            self.group = 'CHEM-L=$'
        else:
            self.group = 'IJ-AVG-$'

        species_in_csv=False
        for row in species_csv:
            if (str(self.name) == row[0].strip()):
                species_in_csv=True
                self.formula   = row[1]
                self.InChI     = row[2]
                self.smiles    = row[3]
                try:
                    self.RMM       = float(row[4])
                except:
                    self.RMM        = 0.0

                if row[5].isspace() or row[5] == "": # Check if the Latex is defined (not whitespace or empty)
                    self.Latex = name
                else: 
                   self.Latex     = row[5]
        if not species_in_csv:               
            print "Warning:"
            print "In MChem_tools.py class Species: Species " + name + " not found in species CSV file"
            return
    def help(self):
        '''
        Another way to get help for the class.
        '''
        help(self)
        return
                                                                    
                                                                                   

