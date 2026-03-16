import TRUST_Post_Loader as tpl

file = "upwind.cgns"
cgns = tpl.TRUST_Post_Loader(file)

times = cgns.getTimes()
print (times)

field_names = cgns.getFieldNames()
print (field_names)

fld = cgns.getFieldDouble("PSCAL_VIT_VIT_SOM_SOM_dom_SOM", -1)
print (fld.getArray())

mesh = cgns.getMesh("dom_ELEM", -1)

# from mc
print(mesh.getNumberOfNodes())
