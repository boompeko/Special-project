#!/usr/bin/env python
# coding: utf-8

# ## Netflix data Analysis

# ## imports

# In[1]:


#importing neccesary library and data 
import numpy as np
import pandas as pd 

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.dates as mdates
get_ipython().run_line_magic('matplotlib', 'inline')

import datetime as dt


# ## Process

# In[2]:


#importing dataset
df= pd.read_csv('netflix1.csv')
df


# In[3]:


#removed all duplicates

# dropped ALL duplicate values
df.drop_duplicates(subset ="show_id",
                     keep = False, inplace = True)
df.head()


# In[4]:


df.info()


# In[5]:


#changing datatype of date_added column to datetime
df["date_added"] = pd.to_datetime(df["date_added"])
df.info()


# ## Analyze

# In[6]:


# added month_added column
df['month_added'] = df['date_added'].dt.month_name()
df


# In[7]:


# added year_added column
df['year_added'] = df['date_added'].dt.year
df


# In[8]:


# added day_added column
df['day_added']=df['date_added'].dt.day_name()
df


# In[9]:


#types of show on netflix
types =df.groupby(['type',])[ 'type'].count().reset_index(name='count')
types=types.set_index('type')
types


# In[10]:


#grouped by directors of show and type
show_director= df.groupby(['director','type'])[ 'director'].count().reset_index(name='show_count')
show_director


# In[11]:


#top 10 directors
top10_directors=show_director.query("`show_count` >= 12")
top10_directors


# In[12]:


#grouped by contry and type
show_origin= df.groupby(['country','type'])[ 'type'].count().reset_index(name='show_count')
show_origin
#


# In[13]:


#write query to find top20 country
top20_country=show_origin.query("`show_count` >= 81")
top20_country


# 

# In[14]:


#grouped by listed_In and type
show_genre= df.groupby(['listed_in','type'])[ 'type'].count().reset_index(name='show_genre_count')
show_genre


# In[15]:


#write query to find top20 genre
top20_genre=show_genre.query("`show_genre_count` >= 110")
top20_genre


# In[16]:


#grouped by rating and type
show_rating= df.groupby(['rating','type'])[ 'rating'].count().reset_index(name='ratings_count')
show_rating=show_rating.set_index('rating')
show_rating


# In[17]:


# filtered  by type = movie

rating_movie= df.groupby(['rating','type'])['rating'].count().reset_index(name='ratings_count')
rating_movie=rating_movie.set_index('rating')
rating_movie

filter = rating_movie["type"]=="Movie"
  
# filtering data
rating_movie.where(filter, inplace = True)
rating_movie= rating_movie.dropna()
rating_movie


# In[18]:


# filtered  by type = TV show

rating_TV= df.groupby(['rating','type'])['rating'].count().reset_index(name='ratings_count')
rating_TV=rating_TV.set_index('rating')
rating_TV

filter = rating_TV["type"]=="TV Show"
  
# filtering data
rating_TV.where(filter, inplace = True)
rating_TV= rating_TV.dropna()
rating_TV


# In[19]:


dsc= df.groupby(['duration','type'])[ 'type'].count().reset_index(name='dsc')
dsc
#SORT
dsc.sort_values(by=['type'], ascending=False)


# In[20]:


#grouped by month_added and type
release_month= df.groupby(['month_added','type'])[ 'type'].count().reset_index(name='release_month')
release_month


# In[21]:


#grouped by year_added and type
release_year= df.groupby(['year_added','type'])[ 'type'].count().reset_index(name='release_count')
release_year


# In[22]:


#grouped by day_added and type
release_Day= df.groupby(['day_added','type'])[ 'type'].count().reset_index(name='release_Day')
release_Day


# ## Share

# In[23]:


# percentage of types of show
colors = ['#ff9999','#66b3ff']
types.plot.pie(y='count',autopct='%.1f%%', shadow=True, legend= 'type' , figsize=(6,6),colors=colors)
plt.title('types_of_show', fontsize=20)


# 

# In[24]:


sns.barplot(x =top20_country.reset_index()['country'], y=top20_country.reset_index()['show_count'],
            hue =  top20_country.reset_index()['type']);
plt.title(' overall_top20_country', fontsize=16)
plt.xlabel('country', fontsize=16);
plt.ylabel('show_count', fontsize=16);
sns.set(rc = {'figure.figsize':(20,10)})
plt.xticks(rotation = 45)


# In[25]:


sns.barplot(x =top20_genre.reset_index()['listed_in'], y=top20_genre.reset_index()['show_genre_count'],
            hue =  top20_genre.reset_index()['type']);
plt.title(' overall_top20_genre', fontsize=16)
plt.xlabel('listed_in', fontsize=16);
plt.ylabel('show_genre_count', fontsize=16);
sns.set(rc = {'figure.figsize':(20,10)})
plt.xticks(rotation = 65)


# In[26]:


explode = (0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,0.05, 0.05, 0.05, 0.05, 0.05,0.05, 0.05,0.05)
show_rating.plot.pie(y='ratings_count',autopct='%.1f%%', figsize=(13,17), pctdistance= 0.80,  explode=explode)
plt.title('overall_ratings distribution', fontsize=20)
centre_circle = plt.Circle((0, 0), 0.70, fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
plt.legend(ncol=15, loc="upper center")


# In[27]:


explode = (0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 )
rating_movie.plot.pie(y='ratings_count',autopct='%.1f%%', figsize=(13,17), pctdistance= 0.80,  explode=explode)
plt.title('Movie_ratings distribution', fontsize=20)
centre_circle = plt.Circle((0, 0), 0.82, fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
plt.legend(ncol=15, loc="upper center")


# In[28]:


explode = (0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,)
rating_TV.plot.pie(y='ratings_count',autopct='%.1f%%', figsize=(13,17), pctdistance= 0.80,  explode=explode)
plt.title('ratings distribution_TV', fontsize=20)
centre_circle = plt.Circle((0, 0), 0.70, fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
plt.legend(ncol=15, loc="upper center")


# In[29]:


sns.barplot(x =release_month.reset_index()['month_added'], y= release_month.reset_index()['release_month'],
            hue =  release_month.reset_index()['type']);
plt.title('release_month' , fontsize=16)
plt.xlabel('month', fontsize=16);
plt.ylabel('show_count', fontsize=16);
sns.set(rc = {'figure.figsize':(12,6)})


# In[30]:


sns.set_style("whitegrid")
# plot boxplot
gfg = sns.lineplot(x ="year_added", y ="release_count", hue="type" ,style="type", markers=True, data = release_year)
 
# add label to the axis and label to the plot
gfg.set(xlabel ="year", ylabel = "release_count", title ='release_per_year')


# In[31]:


sns.barplot(x =release_Day.reset_index()['day_added'], y= release_Day.reset_index()['release_Day'],
            hue =  release_Day.reset_index()['type']);
plt.title('release_day' , fontsize=16)
plt.xlabel('day_of_week', fontsize=16);
plt.ylabel('show_count', fontsize=16);
sns.set(rc = {'figure.figsize':(12,6)})


# ## key findings

# *** A pie chart comparing Netflix movie uploads to TV show uploads from 2008 to 2021 reveals a 39.4% increase in movie uploads.**
# 
# 
# *** Netflix has the most movies from the United States, followed by India in second place and the United Kingdom in third place. It also has the most TV episodes from the United States, followed by Pakistan in second place and the United Kingdom in third place.**
# 
# *** Netflix has the most titles in the "Dramas, International Movies" genre, followed byCrime TV Shows, International TV Shows, TV Dramas' in second place and Stand-Up Comedy in third, according to the overall top 20 genre barchart. and one of the key findings is that the majority of TV programmes are intended for kids tv.**
# 
# *** In the overall rating distribution donut chart, we can see that 36.5% of the shows have TV-MA ratings, indicating that the majority of the shows are for mature audiences, with TV-14 ratings coming in second. Parental Guidelines denotes content, and R rating at the third R classification denotes that the film is not suitable for minors to watch due to violence, offensive language, or sexual activity.**
# 
# 
# *** The movie rating distribution donut chart reveals that 33.7% of movies are rated TV-MA. A TV-MA rating indicates that the show is meant for mature audiences. The second most shows are classified TV-14, which means that Parental Guidelines indicates content for mature audiences, and the third most shows are rated R, which implies that they are not acceptable for children to watch due to violence, foul language, or sexual activity.**
# 
# *** The TV rating distribution donut chart reveals that 42.9% of moviesTV show are rated TV-MA, which means that the broadcasts are meant for mature audiences. The second most shows are rated TV-14, which implies that Parental Guidelines indicates content for mature audiences, and the third most series are rated TV-PG, which denotes under parental guidance.**
# 
# *** The release month barchart indicates that Netflix routinely releases new shows throughout the year.**
# 
# *** In the release per year linechart, we can observe that Netflix began to add shows in large quantities from 2014.**
# 
# *** The majority of the shows on Netflix are released on Friday, as can be seen in the release day barchart.**

# 

# In[ ]:




