# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 23:53:23 2016

@author: Doug
"""

from panoptes_client import SubjectSet, Subject, Project, Panoptes
import glob
#print glob.glob("c:\\Users\\Doug\\LSST\\sub_sets\\Asteroids\\*.jpg")[0][38:-4]

#Connect to Zooniverse and Find Project
Panoptes.connect(username='dougbrn', password='roscoe282306')
project = Project.find(3356)


#Ping Zooniverse to send data. Data sent in an email.    
#project.get_export("classifications",generate=True, wait=True)

#Create new subject set. Name must be unique.   
subject_set = SubjectSet()
subject_set.links.project = project
subject_set.display_name = 'Asteroids (scaled)'
subject_set.save()

#Load an image.
for image in glob.glob("c:\\Users\\Doug\\LSST\\LSST-Zoo\\sub_sets\\Asteroids\\*.jpg"):
    subject = Subject()
    subject.links.project = project
    subject.add_location(image)
    subject.metadata['image_title'] = image[47:-4]
    subject.save()
    subject_set.add(subject)
# You can set whatever metadata you want, or none at all
#subject.metadata['image_id'] = 12345
    

#SubjectSet.add() #can take a list of Subjects, or just one.

