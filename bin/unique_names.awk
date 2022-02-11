#!/usr/bin/awk -f
{file_list[NR]=$1} 
END{
	# For an unlimited loop
	for(i=0; i<=1000; i++){
		# For each filename
		for(j in file_list){
			# Get the substring that is the last i characters
			sub_string = substr(file_list[j], length(file_list[j])-(i));
			# Increase the count for the number of times that has increased
			current_chars[sub_string]++;
		}
		# Once the current_chars has been populated
		if(current_chars[sub_string]==NR){
			# Then we need to keep going onto the next i
			delete current_chars;
		}else{
			# Then the previous i was the length that we want to remove from each name
			for(k in file_list){
				print substr(file_list[k],1, length(file_list[k])-(i));
			};
			break;
		}
	}
}
