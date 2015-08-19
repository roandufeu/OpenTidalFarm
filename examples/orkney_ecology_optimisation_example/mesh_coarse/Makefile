mesh:
	@echo "Generating mesh..."
	#gmsh -2 coast.geo
	dolfin-convert -o xml coast.msh coast.xml
	convert_to_new_xml.py coast.xml

mesh_coarse_domain_farm_2.xml: coast.xml
	python create_mesh_domain_coarse.py coast 2 $@
