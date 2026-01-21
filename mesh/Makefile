# File geo
GEO = mesh_fsi_geo.geo

# File mesh 
MESH = mesh_fsi.msh

# gmsh 
GMSH = gmsh

# Mesh options
# -2  -> genera mesh 2D
OPTIONS = -2

.PHONY: all clean mesh

all: mesh

mesh: $(MESH)

$(MESH): $(GEO)
	$(GMSH) $(OPTIONS) $(GEO) -o $(MESH)

clean:
	rm -f $(MESH)
