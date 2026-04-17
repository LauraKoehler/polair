# polair

## preprocessing

```
polair preprocessing -f <flight number> -c <config file>
```

## noseboom

```
polair noseboom -f <flight number> -c <config file>
```

## tbird

```
polair tbird -f <flight number> -c <config file>
```

# Find corrupted lines in DMS files

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