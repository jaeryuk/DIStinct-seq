import pandas as pd
import numpy as np
import sys

##file path
sam_file = sys.argv[2] + "/" + sys.argv[1] + "_final.sam"

header=('read_ID', 'sam_flag', 'chr', 'position', 'MAPQ', 'CIGAR', 'mate_chr', 'mate_position', 'frag_length', 'base_call', 'base_call_Q', 'tag1', 'tag2', 'tag3', 'tag4', 'tag5', 'tag6', 'tag7', 'tag8')
df = pd.read_csv(sam_file, sep="\t", names=header)

#create position ID column
df["position_ID"] = df["chr"] + ":" + df.apply(lambda x: -x["position"] if (x['frag_length'] < 0)  else x["position"] , axis=1).astype(str)

#create read length coloumn
df["read_length"] = df["CIGAR"].str.extractall('(\d+)').astype(int).groupby(level=0).sum().squeeze()

#create strand column
df['strand'] = df.apply(lambda x: "-" if x['frag_length'] < 0 else "+" , axis=1)

#create alternative alignment position dataframe
df_XA = df[df["tag1"].str.contains("XA")==True]
XA_edit1 = df_XA['tag1'].str.split(":|;|,",expand=True)

df_XA["secondary_position_1"] = XA_edit1[2] + ":" + XA_edit1[3].astype(float).astype('Int64').astype(str)
df_XA["secondary_position_2"] = XA_edit1[6] + ":" + XA_edit1[7].astype(float).astype('Int64').astype(str)
df_XA["secondary_position_3"] = XA_edit1[10] + ":" + XA_edit1[11].astype(float).astype('Int64').astype(str)
df_XA["secondary_position_4"] = XA_edit1[14] + ":" + XA_edit1[15].astype(float).astype('Int64').astype(str)
df_XA["secondary_position_5"] = XA_edit1[18] + ":" + XA_edit1[19].astype(float).astype('Int64').astype(str)

df_XA = df_XA.replace({":<NA>": np.NaN})

df_XA_edit1 = df_XA[['position_ID', 'secondary_position_1', 'secondary_position_2','secondary_position_3','secondary_position_4','secondary_position_5']]

##drop duplicated positioin
df_XA_edit2 = df_XA_edit1.drop_duplicates()
df_XA_edit2.reset_index(inplace=True, drop=True)

#secondary positions corresponding primary position
df_XA_edit3 = pd.melt(df_XA_edit2, id_vars=["position_ID"],
                  var_name="secondary_number", value_name="secondary_position")

df_XA_edit4 = df_XA_edit3.dropna()
df_XA_edit5 = df_XA_edit4.sort_values(by=["position_ID", "secondary_number"])

#list of primary position ID sorted by count
df_XA_edit5_order = df.groupby(['position_ID'])['read_ID'].count().reset_index(name='Count').sort_values(['Count'], ascending=False)["position_ID"].to_list()

#sort by count of position ID
df_XA_edit5['position_ID'] = pd.Categorical(df_XA_edit5['position_ID'], df_XA_edit5_order)
df_XA_edit6 = df_XA_edit5.sort_values('position_ID')
df_XA_edit6.reset_index(inplace=True, drop=True)

#only keep primary positions not in seconday positions of more abundant primary position
secondary_position_collection = []
result =[]
for index, row in df_XA_edit6.iterrows():
    appeared_secondary_position = row['secondary_position']
    secondary_position_collection.append(appeared_secondary_position)

    if row['position_ID'] not in (secondary_position_collection):
        result.append([ row['position_ID'], row['secondary_position']])

#rename column
df_XA_edit7 = pd.DataFrame(result, columns=["position_ID", "secondary_position"])
df_XA_edit7["primary_count"] = df_XA_edit7.groupby("position_ID")['secondary_position'].transform('size')
df_XA_edit7["secondary_count"] = df_XA_edit7.groupby("secondary_position")['position_ID'].transform('size')
df_XA_edit8 = df_XA_edit7.drop_duplicates().reset_index(drop=True)
df_XA_edit9 = df_XA_edit8.rename(columns={'position_ID': 'primary_position_ID', 'secondary_position': 'position_ID'})

#only keep positions of more probable primary positions
df_XA_edit10 = df_XA_edit9[df_XA_edit9["primary_count"] > df_XA_edit9["secondary_count"]]

#only keep secondary positions not in primary positions
primary_position_collection = []
result =[]
for index, row in df_XA_edit10.iterrows():
    appeared_primary_position = row['primary_position_ID']
    primary_position_collection.append(appeared_primary_position)

    if row['position_ID'] not in (primary_position_collection):
        result.append([ row['primary_position_ID'], row['position_ID'], row['primary_count'], row['secondary_count']])

#rename column
df_XA_edit11 = pd.DataFrame(result, columns=["primary_position_ID", "position_ID", "primary_count", "secondary_count"])

##remove duplicated postions
mask = df_XA_edit11['position_ID'].duplicated(keep="first")
df_XA_edit12 = df_XA_edit11[~mask]

##replace secondary postions to correspoding primary postions
df_edit1 = df.join(df_XA_edit12.set_index('position_ID')['primary_position_ID'], on='position_ID')

df_edit1['position2'] = np.where(  df_edit1['strand'] == "-"  , df_edit1['position'].abs() + df_edit1['read_length'] - 2, df_edit1["position"].abs())
df_edit1['position_ID2'] = df_edit1["position_ID"].str.split(":",expand=True)[0] + ":" + df_edit1["position2"].astype(str)

df_edit1["new_position_ID"] = np.where(df_edit1["primary_position_ID"].isnull() , df_edit1["position_ID"], df_edit1["primary_position_ID"])
df_edit1["new_position"] = df_edit1["new_position_ID"].str.split(":", expand=True)[1].astype(int)

#create distinct postion considering strands
df_edit1['distinct_position'] = np.where(  df_edit1['new_position'] < 0  , df_edit1['new_position'].abs() + df_edit1['read_length'] - 2, df_edit1["new_position"].abs())
df_edit1['distinct_position_ID'] = df_edit1["new_position_ID"].str.split(":",expand=True)[0] + ":" + np.where(  df_edit1['new_position'] < 0  , "-", "+") + df_edit1["distinct_position"].astype(str)

#break point extract
break_point1 = df_edit1.value_counts(["distinct_position_ID"]).to_frame('count').reset_index()

break_point1["chr"] = break_point1["distinct_position_ID"].str.split(":", expand=True)[0]
break_point1["position"] = break_point1["distinct_position_ID"].str.split(":", expand=True)[1]

#ambiguous position create
break_point1["ambiguous1"] = break_point1["position"].astype(int) - 3
break_point1["ambiguous2"] = break_point1["position"].astype(int) - 2
break_point1["ambiguous3"] = break_point1["position"].astype(int) - 1
break_point1["ambiguous4"] = break_point1["position"].astype(int) + 1
break_point1["ambiguous5"] = break_point1["position"].astype(int) + 2
break_point1["ambiguous6"] = break_point1["position"].astype(int) + 3

break_point1["ambiguous1"] = break_point1["chr"] + ":" + break_point1["ambiguous1"].astype(str)
break_point1["ambiguous2"] = break_point1["chr"] + ":" + break_point1["ambiguous2"].astype(str)
break_point1["ambiguous3"] = break_point1["chr"] + ":" + break_point1["ambiguous3"].astype(str)
break_point1["ambiguous4"] = break_point1["chr"] + ":" + break_point1["ambiguous4"].astype(str)
break_point1["ambiguous5"] = break_point1["chr"] + ":" + break_point1["ambiguous5"].astype(str)
break_point1["ambiguous6"] = break_point1["chr"] + ":" + break_point1["ambiguous6"].astype(str)

break_point2 = break_point1.drop(["count", "chr", "position"], axis=1)
break_point3 = pd.melt(break_point2, id_vars=["distinct_position_ID"], var_name="ambiguous", value_name="ambiguous_position")

break_point1_order = break_point1["distinct_position_ID"].to_list()
pd.Categorical(break_point3['distinct_position_ID'], break_point1_order)
break_point3["distinct_position_ID"] = pd.Categorical(break_point3['distinct_position_ID'], break_point1_order)
break_point4 = break_point3.sort_values(['distinct_position_ID', 'distinct_position_ID']).reset_index(drop=True)

#only keep distinct_position not in ambiguous_position
ambiguous_position_collection = []
result =[]
for index, row in break_point4.iterrows():
    appeared_ambiguous_position = row['ambiguous_position']
    ambiguous_position_collection.append(appeared_ambiguous_position)

    if row['distinct_position_ID'] not in (ambiguous_position_collection):
        result.append([ row['distinct_position_ID'], row['ambiguous_position']])

break_point5 = pd.DataFrame(result, columns=["distinct_position_ID", "ambiguous_position"])

##remove duplicated ambiguous position
mask = break_point5['ambiguous_position'].duplicated(keep='first')
break_point6 = break_point5[~mask]
break_point6 = break_point6.rename(columns={"distinct_position_ID":"new_distinct_position_ID"})

#merge with sam file
df_edit2 = pd.merge(df_edit1, break_point6, how='left', left_on='distinct_position_ID', right_on='ambiguous_position')
df_edit2["distinct_position_ID2"] = np.where(df_edit2["new_distinct_position_ID"].isnull() , df_edit2["distinct_position_ID"], df_edit2["new_distinct_position_ID"])

#break point count and mean read length
break_point7 = df_edit2.groupby('distinct_position_ID2') \
       .agg({'read_ID':'size', 'read_length':'mean'}) \
       .rename(columns={'read_ID':'count','read_length':'mean_read_length'}) \
       .reset_index()

#only keep reads whth mean read length above 50 to exclude false positive
break_point7 = break_point7[break_point7["mean_read_length"] >= 50].sort_values("count", ascending=False)

#only keep reads with count above 10 to exclude false positive
# break_point8 = break_point7[break_point7["count"] >= 10]

#exclude decoy or alt chromosomes
searchfor = ['_', '\*']
break_point9 = break_point7[~break_point7["distinct_position_ID2"].str.contains('|'.join(searchfor))]


break_point9["final_position_ID"] = break_point9["distinct_position_ID2"].str.split(":", expand=True)[0] + ":" +  break_point9["distinct_position_ID2"].str.split(":", expand=True)[1].astype(int).abs().astype(str)
break_point9["final_strand"] = np.where(break_point9["distinct_position_ID2"].str.split(":", expand=True)[1].astype(int) < 0, "-", "+")

#create final breakpoint file
break_point9["chr"] = break_point9["final_position_ID"].str.split(":", expand=True)[0]
break_point9["start"] =break_point9["final_position_ID"].str.split(":", expand=True)[1]
break_point9["end"] = break_point9["final_position_ID"].str.split(":", expand=True)[1].astype(int) + 1
break_point9["null"] = "."

break_point10 = break_point9[["chr", "start", "end", "count", "null", "final_strand"]]

#save
break_point10.to_csv(sys.argv[2] + "/" + sys.argv[1] + "_breakpoint.txt", sep="\t", header=False, index=False)
