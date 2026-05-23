# File geo
GEO = mesh_fsi_geo.geo

# File mesh 
MESH = mesh_fsi.msh

# gmsh 
GMSH = gmsh

# Mesh options
# -2        -> generate a 2D mesh
# -format   -> enforce MSH v2 for deal.II GridIn::read_msh
# -nopopup  -> avoid interactive GUI mode
OPTIONS = -2 -format msh2 -nopopup

.PHONY: all clean mesh

all: mesh

mesh: $(MESH)

$(MESH): $(GEO)
	$(GMSH) $(OPTIONS) $(GEO) -o $(MESH)

clean:
	rm -f $(MESH)
