#!/bin/rc
awk '
NR == 1{
	n = NF
	for(i = 1; i <= n; i++)
		max[i] = min[i] = $i
	next
}
{
	for(i = 1; i <= n; i++)
		if($i > max[i])
			max[i] = $i
		else if($i < min[i])
			min[i] = $i
}
END{
	for(i = 1; i <= n; i++)
		print min[i], max[i]
}
' $*
