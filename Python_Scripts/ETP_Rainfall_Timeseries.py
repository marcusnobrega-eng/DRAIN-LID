from datetime import datetime
from meteostat import Point, Daily
import numpy as np
import pandas as pd
import os

# === Correct Hargreaves ETP function ===
def hargreaves_et0(tavg, tmax, tmin, lat_deg, doy):
    phi = np.deg2rad(lat_deg)
    Gsc = 1360  # Solar constant [J/m2/s]
    Lambda = 2.45  # Latent heat [MJ/kg]

    delta = 0.409 * np.sin((2 * np.pi * doy / 365) - 1.39)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)
    omega_s = np.arccos(-np.tan(phi) * np.tan(delta))

    Ra = (1 / np.pi) * Gsc * dr * (
        omega_s * np.sin(phi) * np.sin(delta) +
        np.cos(phi) * np.cos(delta) * np.sin(omega_s)
    )

    Ra_mm_day = (Ra * 86400) / (Lambda * 1e6)
    Tdiff = np.maximum(tmax - tmin, 0)
    et0 = 0.0023 * (tavg + 17.8) * np.sqrt(Tdiff) * Ra_mm_day
    return et0

# === Time period ===
start = datetime(1995, 1, 1)
end = datetime(2025, 1, 1)

# === City coordinates ===
airports = {
    'New York City': Point(40.6999, -74.1755, 4),
    'Miami': Point(25.90647, -80.27489, 3),
    'Phoenix': Point(33.43251, -112.0116, 346),
    'San Francisco': Point(37.6203, -122.38435, 4)
}
latitudes = {
    'New York City': 40.6999,
    'Miami': 25.90647,
    'Phoenix': 33.43251,
    'San Francisco': 37.6203
}

# === Create output directory ===
os.makedirs("ETP_Hargreaves_Exports", exist_ok=True)

# === Summary stats storage ===
summary_stats = []

# === Loop through each city ===
for city, location in airports.items():
    print(f"\n📡 Fetching data for: {city}")
    df = Daily(location, start, end).fetch()
    df = df[['tmin', 'tmax', 'prcp']].copy()
    df.reset_index(inplace=True)
    df.dropna(inplace=True)

    # === Compute ETP ===
    df['Julian'] = df['time'].dt.dayofyear
    df['Tavg'] = (df['tmax'] + df['tmin']) / 2
    lat = latitudes[city]
    df['ETP_mm_day'] = hargreaves_et0(df['Tavg'], df['tmax'], df['tmin'], lat, df['Julian'])

    # === Convert to mm/h ===
    df['Rainfall_mm_h'] = df['prcp'] / 24
    df['ETP_mm_h'] = df['ETP_mm_day'] / 24
    df['Inflow_Hydrograph_mm_h'] = df['ETP_mm_h'] - df['Rainfall_mm_h']

    # === Time in minutes since start ===
    df['Time (min)'] = ((df['time'] - df['time'].iloc[0]).dt.total_seconds() / 60).astype(int)

    # === Save formatted time series ===
    output_df = df[['Time (min)', 'Rainfall_mm_h', 'ETP_mm_h', 'Inflow_Hydrograph_mm_h']]
    output_df.columns = ['Time (min)', 'Rainfall (mm/h)', 'ETP (mm/h)', 'Inflow Hydrograph (mm/h)']
    filename = f"ETP_Hargreaves_Exports/{city.replace(' ', '_')}_ETP_Timeseries.csv"
    output_df.to_csv(filename, index=False)
    print(f"✅ Exported to: {filename}")

    # === Annual statistics ===
    df['Year'] = df['time'].dt.year
    annual_rain = df.groupby('Year')['prcp'].sum()
    annual_etp = df.groupby('Year')['ETP_mm_day'].sum()

    mean_rain = annual_rain.mean()
    std_rain = annual_rain.std()
    mean_etp = annual_etp.mean()
    std_etp = annual_etp.std()

    aridity_index = mean_rain / mean_etp if mean_etp > 0 else np.nan

    summary_stats.append({
        'City': city,
        'Mean Annual Rainfall (mm)': f"{mean_rain:.0f} ± {std_rain:.0f}",
        'Mean Annual ETP (mm)': f"{mean_etp:.0f} ± {std_etp:.0f}",
        'Aridity Index (Rainfall / ETP)': round(aridity_index, 3)
    })

# === Export summary CSV ===
summary_df = pd.DataFrame(summary_stats)
summary_path = "ETP_Hargreaves_Exports/Climate_Summary_Statistics.csv"
summary_df.to_csv(summary_path, index=False, encoding='utf-8-sig')
print(f"\n📊 Climate summary statistics saved to: {summary_path}")
