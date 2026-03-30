import math

### add percentage dispersion term (such as /3 *0.3 for 1 value and keeping 1 /3 and then allowing one to go higher)
### add correct key ordering
### calculate standard deviation
### calculate rank-based standard deviation

log_values = [0.2, 0.3, 0.5, 0.8, 1, 1.2, 2.4, 2, 3.5, 4, 5, 6, 8, 12, 15, 16, 19]
log_transform = lambda x: 5*math.log(1+x, 10)
reverse_log = lambda x: (10**(x/5))-1
original_values = [reverse_log(x)for x in log_values]
# print(original_values)
# print([round(log_transform(x),2) for x in original_values])
max_values_log = max(log_values) + 1
max_values_original = max(original_values) + 1

spread_method = 1 # 1
only_higher = True

methods = [
# "linear_square",
"linear",
"fractional",
# "fractional_square"
]
thresholds = ["none", "local"]

to_compare_with = [1,2,3]
values_to_check = [
# [1],
# [1,2],
# [1,5],
# [4],
# [6],
# [6,8],
# [2,10,11],
[10,11,15],
[12],
[12,13,14],
[14,15],
[4,5,6],
[15],
]
values_ = [log_values[x] for positions in values_to_check for x in positions]

log_values_new = log_values.copy()
original_values_new = original_values.copy()
log_dict_ranks = {}
original_dict_ranks = {}
log_dict_values = {}
original_dict_values ={}

for method in methods:
	log_dict_ranks[method] = {}
	original_dict_ranks[method] = {}
	log_dict_values[method] = {}
	original_dict_values[method] = {}

for compare_with_value in to_compare_with:
	for positions in values_to_check:
		for threshold in thresholds:
			if threshold == "local":
				values = log_values
				max_values = max_values_log
			elif threshold == "none":
				values = original_values	
				max_values = max_values_original	
			for method in methods:
				values_ = [values[x] for x in positions]
				if method == "linear":
					left_hand_side = [max_values - x for x in values_]
				elif method == "linear_square":
					left_hand_side = [(max_values - x)**2 for x in values_]
				elif method == "fractional":
					left_hand_side = [1/x for x in values_]
				elif method == "fractional_square":
					left_hand_side = [1/(x**2) for x in values_]
				if len(left_hand_side) == compare_with_value:
					continue
				elif only_higher and len(left_hand_side) > compare_with_value:
					continue
				print(method, threshold, values_)
				sum_lhs = sum(left_hand_side)

				if compare_with_value == 1:
					right_hand_side = [sum_lhs]
				else:
					if spread_method != 1:
						temp_right_hand_side_values = []
						initial_value = sum_lhs/compare_with_value
						for idx in range(1, compare_with_value):
							spread_to_use = 1 - (spread_method * idx)
							if spread_to_use < 0:
								spread_to_use = max(spread_method * idx, compare_with_value - 0.5)
							temp_right_hand_side_values.append(initial_value * spread_to_use)
						sum_temp_rhs = sum(temp_right_hand_side_values)
						difference = sum_lhs - sum_temp_rhs
						temp_right_hand_side_values.append(difference)
						right_hand_side = temp_right_hand_side_values.copy()				

					else:
						right_hand_side = [sum_lhs/compare_with_value for idx in range(compare_with_value)]
				if method == "linear":
					right_hand_side_numbers = [max_values - rhs_value for rhs_value in right_hand_side]
					if max(x for x in right_hand_side) > max_values:
						continue
				elif method == "linear_square":
					if max(math.sqrt(x) for x in right_hand_side) > max_values:
						continue
					right_hand_side_numbers = [max_values - math.sqrt(rhs_value) for rhs_value in right_hand_side]
				elif method == "fractional":
					right_hand_side_numbers = [1/x for x in right_hand_side]
				elif method == "fractional_square":
					right_hand_side_numbers = [1/math.sqrt(x) for x in right_hand_side]
				if threshold == "local":
					log_values_new.extend(right_hand_side_numbers)
					original_values_new.extend(round(reverse_log(x),3) for x in right_hand_side_numbers)
				elif threshold == "none":
					# print(compare_with_value)
					log_values_new.extend(log_transform(x) for x in right_hand_side_numbers)
					# print(log_values_new)
					original_values_new.extend(right_hand_side_numbers)

log_values_new = sorted(list(set(log_values_new)))
original_values_new = sorted(list(set(original_values_new)))

# log_values_new.sort()
# original_values_new.sort()
for compare_with_value in to_compare_with:
	for positions in values_to_check:
		for threshold in thresholds:
			if threshold == "local":
				values = log_values
				max_values = max_values_log
			elif threshold == "none":
				values = original_values	
				max_values = max_values_original	
			for method in methods:
				values_ = [values[x] for x in positions]
				if method == "linear":
					left_hand_side = [max_values - x for x in values_]
				elif method == "linear_square":
					# print(values_,2)
					left_hand_side = [(max_values - x)**2 for x in values_]
					# print(left_hand_side,3)
				elif method == "fractional":
					left_hand_side = [1/x for x in values_]
				elif method == "fractional_square":
					left_hand_side = [1/(x**2) for x in values_]
				if len(left_hand_side) == compare_with_value:
					continue
				elif only_higher and len(left_hand_side) > compare_with_value:
					continue
				# print(left_hand_side)
				sum_lhs = sum(left_hand_side)
				if compare_with_value == 1:
					right_hand_side = [sum_lhs]
				else:
					if spread_method != 1:
						temp_right_hand_side_values = []
						initial_value = sum_lhs/compare_with_value
						for idx in range(1, compare_with_value):
							spread_to_use = 1 - (spread_method * idx)
							if spread_to_use < 0:
								spread_to_use = max(spread_method * idx, compare_with_value - 0.5)
							temp_right_hand_side_values.append(initial_value * spread_to_use)
						sum_temp_rhs = sum(temp_right_hand_side_values)
						difference = sum_lhs - sum_temp_rhs
						temp_right_hand_side_values.append(difference)
						right_hand_side = temp_right_hand_side_values.copy()				

					else:
						right_hand_side = [sum_lhs/compare_with_value for idx in range(compare_with_value)]
				if method == "linear":
					right_hand_side_numbers = [max_values - rhs_value for rhs_value in right_hand_side]
					if max(x for x in right_hand_side) > max_values:
						continue
				elif method == "linear_square":
					if max(math.sqrt(x) for x in right_hand_side) > max_values:
						# print(right_hand_side, 999)
						continue
					right_hand_side_numbers = [max_values - math.sqrt(rhs_value) for rhs_value in right_hand_side]
				elif method == "fractional":
					right_hand_side_numbers = [1/x for x in right_hand_side]
				elif method == "fractional_square":
					right_hand_side_numbers = [1/math.sqrt(x) for x in right_hand_side]
				if threshold == "local":
					left_hand_side_indices = [log_values_new.index(x) for x in values_]
					log_dict_ranks[method][str(left_hand_side_indices)+ "to " + str(compare_with_value)] = [log_values_new.index(x) for x in right_hand_side_numbers]
					# original_dict_ranks[method][str(left_hand_side_indices)] = [original_values_new.index(round(reverse_log(x),3)) for x in right_hand_side_numbers]
					log_dict_values[method][str(values_)+ "to " + str(compare_with_value)] = [x for x in right_hand_side_numbers]
					# original_dict_values[method][str(values_)] = [round(reverse_log(x),3) for x in right_hand_side_numbers]
				elif threshold == "none":
					left_hand_side_indices = [original_values_new.index(x) for x in values_]
					# log_dict_ranks[method][str(left_hand_side_indices)] = [log_values_new.index(log_transform(x)) for x in right_hand_side_numbers]
					original_dict_ranks[method][str(left_hand_side_indices)+ "to " + str(compare_with_value)] = [original_values_new.index(x) for x in right_hand_side_numbers]
					# log_dict_values[method][str(values_)] = [round(log_transform(x),3) for x in right_hand_side_numbers]
					original_dict_values[method][str(values_)+ "to " + str(compare_with_value)] = [x for x in right_hand_side_numbers]
print(log_values)
# print(log_values_new)			
# print(log_dict_ranks)
# print(log_dict_values)
# for key,dict_ in log_dict_ranks.items():
# 	print(key, dict_)
import pandas as pd
df = pd.DataFrame()
for key_, dict_ in log_dict_values.items():
	keys = list(dict_.keys())
	keys.sort()
	print(keys)
	sorted_dict = {i: dict_[i] for i in keys}
	print("\n",key_,sorted_dict)
	for key, values in sorted_dict.items():
		for idx, value in enumerate(values):
			df = df.append({'method': key_, 'key': key, 'values': value, "special_idx": idx}, ignore_index=True)

df = df.explode('values')
print(set(df["method"]))
import plotly.express as px


# Create a new column for displaced category by appending the special_idx to the key
df['displaced_key'] = df['key'].astype(str) + "_" + df['special_idx'].astype(str)

# Convert key to categorical to maintain proper ordering
# df['key'] = pd.Categorical(df['key'], ordered=True)

# Create the scatter plot using the displaced_key for the y-axis
fig = px.scatter(df, x='values', y='key', color='method')

# Update layout to display the original key on the y-axis for better readability
fig.update_layout(
    yaxis_title="Keys (Displaced)",
    title="Scatter Plot with Displaced Y Axis Based on special_idx"
)

fig.show()
# # Create the bar plot
# fig = px.scatter(df, x='values', y='key', color='method', orientation='h')

# # Show the plot
# fig.show()
# print()