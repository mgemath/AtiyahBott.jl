using Test
using AtiyahBott

# first test
println( "Testing O1" )
P = (g,c,w,s,m) -> O1(g,c,w,s,m)^2
res = AtiyahBottFormula(2,2,6-1,P,true,true,false)
if res[1] == 1//1 
    @test true
else
    @test false
end

#second test
println( "Testing Incidency and Contact" )
d = 4
P = Vector{Any}(undef, d+1)
for a in 0:d
    P[a+1] = (g,c,w,s,m) -> Incidency(g,c,w,s,3)^a*
        Incidency(g,c,w,s,2)^(2*(d-a)+1)*
        Contact(g,c,w,s)
end
res = AtiyahBottFormula(3,d,0,P,true,true,false);

if res == [ 1089024, 96512, 9408, 1024, 128 ] 
    @test true
else
    @test false
end

#third test
println( "Testing Hypersurface" )
d = 4 #for other values of d, change this line
P = (g,c,w,s,m) -> Hypersurface(g,c,w,s,5)
res = AtiyahBottFormula(4,d,0,P,true,true,false);

if res == [ 15517926796875//64 ]
    @test true 
else  
    @test false
end

#fourth test
println( "Testing R1" )
d = 4 #for other values of d, change this line
P = (g,c,w,s,m) -> R1(g,c,w,s,1)^2
res = AtiyahBottFormula(1,d,0,P,true,true,false)

if res == [ 1//32 ]
    @test true
else
    @test false
end