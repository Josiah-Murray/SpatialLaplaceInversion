#This file generates the code for the 'CohenSuite' function 
# contained in "TestFunctions.jl"

with open("CohenTestFunctions.txt", "x") as TextFile:


    TextFile.write("""
function CohenSuite( inversionScheme, t)
    exact = zeros(35)
    approx = zeros(35)
""")
    for i in range(1,36):

        TextFile.write(f"""
    exact[{i}] = Cohen{i}_exact(t)
    approx[{i}] = inversionScheme(Cohen{i}, t)
            """)



    TextFile.write("""
    return [exact approx]     
end
        """)