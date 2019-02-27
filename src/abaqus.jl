using Formatting

function writeAbq(filename::String, mesh::Mesh)
    # open file
    open(filename, "w") do afile
        # node header
        write(afile, "*Node\n")
        
        # node coordinates
        for (nodeid, nodecoords) in enumerate(mesh.nodes)
            write(afile, format("{}, {}, {}\n", nodeid, nodecoords...))
        end
        
        # element header
        write(afile, "*Element, type=CPS4R\n")
        
        # elements
        for (elid, element) in enumerate(mesh.elements)
            write(afile, format("{}, {}, {}, {}, {}\n", elid, element...))
        end
    end
end