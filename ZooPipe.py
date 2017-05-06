# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 23:53:23 2016

@author: Doug
"""

from panoptes_client import SubjectSet, Subject, Project, Panoptes
import glob
import csv
#print glob.glob("c:\\Users\\Doug\\LSST\\sub_sets\\Asteroids\\*.jpg")[0][38:-4]

#Connect to Zooniverse and Find Project
Panoptes.connect(username='dougbrn', password='roscoe282306')
project = Project.find(3356)

#Create new subject set. Name must be unique.
subject_set = SubjectSet()
subject_set.links.project = project
subject_set.display_name = 'Asteroids (4/19/17)'
subject_set.save()

#Load an image.
print len(glob.glob("/home/doug/lsst-mount/Zooniverse/cutouts/*.png"))
for image in glob.glob("/home/doug/lsst-mount/Zooniverse/cutouts/*.png"):
    subject = Subject()
    subject.links.project = project
    subject.add_location(image)
    subject.metadata['image_title'] = image
    subject.save()
    subject_set.add(subject)
# You can set whatever metadata you want, or none at all
#subject.metadata['image_id'] = 12345


#SubjectSet.add() #can take a list of Subjects, or just one.
"""
#Ping Zooniverse to send data. Data sent in an email.
export_file = project.get_export("classifications",generate=True, wait=True)
r  = csv.reader(export_file)
for row in r:
    print row[0]
"""
