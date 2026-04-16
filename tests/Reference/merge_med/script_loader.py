import TRUST_Post_Loader as tpl

file = "par_cgns.cgns"
cgns = tpl.TRUST_Post_Loader(file)

times = cgns.getTimes()
print (times)

field_names = cgns.getFieldNames()
print (field_names)

fld = cgns.getFieldDouble("PRESSION_SOM_dom_SOM", -1)
print (fld.getArray())

mesh = cgns.getMesh("dom_ELEM", -1)

# from mc
print(mesh.getNumberOfNodes())
