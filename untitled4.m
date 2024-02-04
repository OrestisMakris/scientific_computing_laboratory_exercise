clear 
a = randn(3, 3, 3); % Example multidimensional tensor a
b = randn(3, 3, 3); % Example multidimensional tensor b
a 

c = ttt_myid(a, b) % Element-wise product (outer product)
cc = ttt(tensor(a),tensor(b))
inner_product = ttt_myid(a, b, 'all') % Inner product

ttt(tensor(a),tensor(b), 1:3 )

size(c)
size(cc)
