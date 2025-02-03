import matplotlib.pyplot as plt

# Read the output.txt file
with open("output.txt", "r") as file:
    lines = file.readlines()

# Extract processor counts and tick counts
processor_counts = []
tick_counts = []

for line in lines:
    if "ticks: " in line:
        ticks = float(line.split(" ")[1].strip())
        tick_counts.append(ticks)
    elif "using " in line:
        processes = int(line.split("using")[1].split()[0])
        processor_counts.append(processes)

# Print extracted data for debugging
print("Processor Counts:", processor_counts)
print("Tick Counts:", tick_counts)

# Create the bar chart
plt.figure(figsize=(8, 6))
plt.bar(processor_counts, tick_counts, color='skyblue')

# Add labels and title
plt.xlabel('Number of Processes')
plt.ylabel('Ticks')
plt.title('Comparison of Ticks for Different Processor Counts')

# Add tick counts on top of each bar
for i, tick in enumerate(tick_counts):
    plt.text(processor_counts[i], tick + 0.5, f'{tick:.2f}', ha='center')

# Show the plot
plt.show()
