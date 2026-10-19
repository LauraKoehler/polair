# How to ...

## Find corrupted lines in DMS files

When downloading the data from the DMS system, it happens that there errors in the files such as broken lines, special charachters, missing entries,... When this is the case, the polair package fails in processing these data. Thus, the corrupted lines have to be removed. Here, we provide some code examples which can be used to quickly identify the corrupted line. In most of the cases, this works:

```python
fn = <corrupted_file>

with open(fn, "rb") as f:
    for i, line in enumerate(f, start=1):
        try:
            line.decode("utf-8")
        except UnicodeDecodeError as e:
            print(f"Unicode error in line {i}: {e}")
            print(line[:200])  # show beginning of bad line
            break
```
This gives the corrupted line number. However, empty lines are not counted. So, if every second line is empty, you need to double the result. If this is not working, try

```python
fn = <corrupted_file>

df = pd.read_csv(fn, header  = 4, sep = r'\s+', names = ["date", "time", f"var"])
bad_idx = [i for i, x in enumerate(df["time"]) if len(str(x)) != 15]
```
bad_idx gives lines with broken timestamps. Duplicated timestamps could also cause problems:

```python
dup_times = df.loc[df["time"].duplicated(), "time"]
```
After removing the corrupted lines manually, the polair commands should work.

## Convert flight segments

The flight segments have to be provided as yaml file to polair. For the campaign documentation on GitHub, it is sometimes useful to have the html code.

### From csv to yaml

Here is the code for converting an csv table to the flight segment yaml file.

```python
fn_in = "FlightSegment.csv" # csv file with columns Flight; Segment; Time; Type; Comment. Use ; for column seperation
df = pd.read_csv(fn_in, sep = ";")
byflight = df.groupby(by = "Flight")

flight_dict = {}
flights = [1,2,3,4,5,6,7,8,9,10,] # flights during the campaign, adjust
for f in flights:
    selected_flight = byflight.get_group(f)
    flight_list = []
    for i in np.arange(len(selected_flight)):
        seg_info = selected_flight.iloc[i]
        if type(seg_info["Comment"]) == type("str"):
            comment = seg_info["Comment"]
        else:
            comment = ""
        segs = {"segment" : int(seg_info["Segment"]), "time": seg_info["Time"], "type": seg_info["Type"], "comment": comment}
        flight_list.append(segs)
    flight_dict[f] = {"segments": flight_list}

flight_segments = {"flights": flight_dict}

fn_out = "flight_segements.yaml" # adjust path
with open(fn_out, 'w') as file:
    yaml.safe_dump(flight_segments, file, sort_keys=False)
```

### From yaml to html

Converting the yaml dictionary into an html table is mainly useful for the campaign GitHub repos. There, we show a table with all flight segments.

```python
import yaml

def generate_html_table(data, output_file):
    html = []

    html.append('<table style="font-size:90%; border-collapse:collapse;">')

    html.append("""
<thead>
<tr>
<th>Flight</th>
<th>Segment</th>
<th>Time</th>
<th>Type</th>
<th>Comment</th>
</tr>
</thead>
<tbody>
""")

    for flight_id, flight_data in data["flights"].items():

        segments = flight_data["segments"]
        rowspan = len(segments)

        for i, seg in enumerate(segments):

            html.append("<tr>")

            # Only write flight number once with rowspan
            if i == 0:
                html.append(f'<td rowspan="{rowspan}">{flight_id}</td>')

            html.append(f"<td>{seg['segment']}</td>")
            html.append(f"<td>{seg['time']}</td>")
            html.append(f"<td>{seg['type']}</td>")
            html.append(f"<td>{seg.get('comment', '')}</td>")

            html.append("</tr>")

    html.append("</tbody></table>")

    with open(output_file, "w") as f:
        f.write("\n".join(html))

    print(f"HTML table written to {output_file}")

fn_flightsegs = "flight_segements.yaml" # update
fn_html_out = "flight_segements.html" # update

with open(fn_flightsegs) as f:
    segments = yaml.safe_load(f)

generate_html_table(segments, fn_html_out)
```