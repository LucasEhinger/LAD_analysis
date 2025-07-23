import csv

MATCH_COLUMNS = ["Target", "Beam Current (uA)", "SHMS angle [deg]", "HMS angle [deg]"]
TIME_COLUMN = "Total Time"

def read_csv(filename):
	with open(filename, newline='') as csvfile:
		reader = csv.DictReader(csvfile)
		rows = list(reader)
	return rows

def make_key(row):
	return (
		row["Target"],
		row["Beam Current (uA)"],
		float(row["SHMS angle [deg]"]),
		float(row["HMS angle [deg]"])
	)

def find_matching_key(accounting_table, key):
	for existing_key in accounting_table:
		try:
			if (
				existing_key[0] == key[0] and
				existing_key[1] == key[1] and
				abs(existing_key[2] - key[2]) <= 1 and
				abs(existing_key[3] - key[3]) <= 1
			):
				return existing_key
		except:
			if all(existing_key[i] == key[i] for i in range(len(key))):
				return existing_key
	return None

def main():
	input_csv = "LAD_runlist.csv"
	rows = read_csv(input_csv)
	accounting_table = {}

	for row in rows:
		key = make_key(row)
		if ":" in row[TIME_COLUMN]:
			h, m, s = map(int, row[TIME_COLUMN].split(":"))
		else:
			h, m, s = 0, 0, 0
		time = (h * 3600 + m * 60 + s) / 3600  # in hours

		match_key = find_matching_key(accounting_table, key)
		if match_key:
			accounting_table[match_key] += time
		else:
			accounting_table[key] = time

	output_csv = "LAD_accounting_output.csv"
	with open(output_csv, "w", newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(MATCH_COLUMNS + [TIME_COLUMN, "Total Charge [mC]"])
		for key, total_time in accounting_table.items():
			try:
				beam_current = float(key[1])
			except ValueError:
				beam_current = 0.0
			total_charge = beam_current * total_time * 3600  # uA * hr * 3600 = uA*s = mC
			writer.writerow([key[0], key[1], key[2], key[3], total_time, total_charge])
	print(f"Accounting table saved to {output_csv}")

if __name__ == "__main__":
	main()