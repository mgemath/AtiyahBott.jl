using Downloads

const current_dir = @__DIR__
const number_trees = [1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]

list_g = open(dirname(current_dir)*"/src/list_trees.txt") 
for n in 3:3#1:4 #d=4

    Dim_dir::String = dirname(current_dir)*"/Data/Dimension$n" #path of the folder containing the colorations
    mkpath(Dim_dir) #create the folder

    seekstart(list_g)
    #open the file containing the list of Prufer sequences of graphs

    for v in 2:(4+1) #run the computation among all graphs with fixed number of vertices
                
        n_trees_nv = number_trees[v - 1]  #we known how many graphs there are with fixed number of vertices
        
        for _ in 1:n_trees_nv  #run the computation for a fixed graph

            str = readline(list_g) #read a new line, we expect a Prufer seq plus the number of automorphisms

            name_file = string(split(str, ',')[1],"0.gz")
            if !(name_file in readdir(Dim_dir))
                url = "https://raw.githubusercontent.com/mgemath/Colorations/main/Dimension$n/$name_file"
                dest = Dim_dir*"/$name_file"
                try
                    Downloads.download(url, dest)
                catch e
                    printstyled(stderr,"ERROR: ", bold=true, color=:red)
                    printstyled(stderr,sprint(showerror,e), color=:light_red)
                    println(stderr)
                    final_state = false
                    break  #end for n_g
                end
            end
            
        end

    end
end
close(list_g)