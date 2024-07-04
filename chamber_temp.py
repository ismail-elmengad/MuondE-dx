import matplotlib.pyplot as plt

import datetime

import datetime

def convert_unix_time_to_minutes(unix_time):
    return unix_time / 60

time = []
temp = []




with open('MDT_temperature_measurements/mdt.bil1a05_ts04.txt', 'r') as file:
    first_line = file.readline().strip()
    first_columns = first_line.split()
    offset = convert_unix_time_to_minutes(float(first_columns[0]))
    for line in file:
        columns = line.split()
        t_time = convert_unix_time_to_minutes(float(columns[0]))
        time.append(t_time - offset)
        temp.append(float(columns[1]))

plt.plot(time, temp, marker='o', color= 'g')
plt.xlabel('Time (minutes)')
plt.ylabel('Temperature (°C)')
plt.title('Temperature during Muon Test Run in BIL1a05 ts00? (June 2023)')

plt.savefig('Temperature BIL1a05ts00?.png')
plt.show()