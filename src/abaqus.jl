using Formatting

function writeabq(filename::String, mesh::Mesh)
    # open file
    open(filename, "w") do afile
        # node header
        write(afile, "*Node\n")
        
        # node coordinates
        for (nodeid, nodecoords) in enumerate(mesh.nodes)
            write(afile, format("{}, {}, {}\n", nodeid, nodecoords...))
        end
        
        # elements
        for (elem_type, element_ids) in mesh.elementsets
            # element header
            write(afile, "*Element, type=$elem_type\n")
            for element_id in element_ids
                element = mesh.elements[element_id]
                write(afile, format("{}, {}, {}, {}, {}\n", element_id, element...))
            end
        end
    end
end