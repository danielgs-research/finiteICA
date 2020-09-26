
a_222_tensor = [ 1 2 3 4 5 6 7 8]
a_222_tensor = reshape(a_222_tensor,[2,2,2])

disp('funny thing that the 2nd element is not the 2nd in first line. ')
disp('but the 2nd row in 1st column')
disp(a_222_tensor(2))
disp('which is different from python.')
disp('we should do a transpose first to get the same behavior.')

disp("but there\'s no transpose for 3 dimensional matrices...")

disp('also, the reshape itself does in that order for us, ')
disp('so if we rely on reshaping and getting the order we want,')
disp(' there would be no problem.')

disp('so...lets being toying with sub2ind and ind2sub')

disp('ind2sub transforms a LINEAR INDEX to a subscript:')


tensor_size = size(a_222_tensor)
a_sub = ind2sub(tensor_size,2)

disp('funny, because if the output was a 3 dimensional list, we would get different')

[d1 d2 d3] = ind2sub(tensor_size,2)

disp('the tensor with d1,d2,d3 as their dimensions')
a_222_tensor(d1,d2,d3)

disp('now, lets use this subscript on sub2ind, to get back to the same_place as')
disp('we where before')

ind = sub2ind(tensor_size,d1,d2,d3)

disp('nice, we got 2. but what if we gave a list to sub2ind?')

ind = sub2ind(tensor_size,[d1,d2,d3])
[d1,d2,d3]






