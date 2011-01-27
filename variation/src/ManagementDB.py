from sqlalchemy.engine.url import URL
from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text, Boolean, Float, Enum
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany
from elixir import setup_all, session, metadata, entities
from elixir.events import after_update
from elixir.events import after_insert
from elixir.events import before_update
from elixir.events import before_insert
from elixir.options import using_table_options_handler  #using_table_options() can only work inside Entity-inherited class.
from sqlalchemy import UniqueConstraint

from datetime import datetime
from datetime import timedelta

from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.orm import scoped_session, sessionmaker

from pymodule.db import ElixirDB
import smtplib
import string
import random
import crypt
import time
import hashlib
import os
#import datetime


__session__ = scoped_session(sessionmaker(autoflush=True,autocommit=True))
#__metadata__ = ThreadLocalMetaData() #2008-11-04 not good for pylon
__metadata__ = MetaData()

class ManagementDB(ElixirDB):
        __doc__ = __doc__
        option_default_dict = ElixirDB.option_default_dict.copy()
        option_default_dict[('drivername', 1,)][0] = 'mysql'
        option_default_dict[('database', 1,)][0] = 'management'
        def __init__(self, **keywords):
                from pymodule.ProcessOptions import ProcessOptions
                ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
                self.setup_engine(metadata=__metadata__, session=__session__, entities=entities)
                """
                from pymodule import ProcessOptions
                ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)

                if getattr(self, 'schema', None):       #for postgres
                        for entity in entities:
                                using_table_options_handler(entity, schema=self.schema)

                metadata.bind = self._url
                setup_all(create_tables=True)   #create_tables=True causes setup_all to call elixir.create_all(), which in turn calls metadata.create_all()
        """




class sys_groups(Entity):
    gid = Field(Integer, primary_key=True)
    groupname = Field(String(30), required=True)
    group_password = Field(String(64),required=True)
    users = OneToMany('sys_users')
    groupusers = ManyToMany('sys_users', tablename='sys_groupmembers')
    using_options(tablename='sys_groups')
    using_table_options(mysql_engine='InnoDB')

class sys_users(Entity):
    uid = Field(Integer, primary_key=True)
    username = Field(String(50))
    realname = Field(String(32))
    shell = Field(String(20),default='/bin/bash')
    password = Field(String(40))
    homedir = Field(String(32),default='/home/data-exchange-home')
    lastchange = Field(String(50))
    min = Field(Integer)
    max = Field(Integer)
    warn = Field(Integer)
    inact = Field(Integer)
    _expire = Field(Integer,colname='expire',synonym='expire')
    confirmed = Field(Boolean,default=False)
    created = Field(DateTime, default=datetime.now)
    email = Field(String(100))
    status = Field(String(1),default = 'N')
    organisation = Field(String(100))
    group = ManyToOne('sys_groups',colname='gid')
    known_partners = ManyToOne('known_partners')
    groups = ManyToMany('sys_groups', tablename='sys_groupmembers')
    using_options(tablename='sys_users')
    using_table_options(mysql_engine='InnoDB')
    
    ADDITIONAL_HOURS = 48
    ADMIN_MAIL = 'uemit.seren@gmail.com'

    def validate_password(self, password):
        """
        Check the password against existing credentials.
        
        :param password: the password that was provided by the user to
            try and authenticate. This is the clear text version that we will
            need to match against the hashed one in the database.
        :type password: unicode object.
        :return: Whether the password is valid.
        :rtype: bool
        
        """
        salt = self.password[:2]
        
        hashed_pass = crypt.crypt(password, salt)
        return self.password == hashed_pass

    def _set_expire(self,expire):
        self._expire = time.mktime(expire.timetuple())
    def _get_expire(self):
        return datetime.fromtimestamp(self._expire)

    expire = property(_get_expire,_set_expire)

    def _set_password_plain(self, password):
        self._password_plain = password
        chars = string.letters + string.digits
        salt = random.choice(chars) + random.choice(chars)
        self.password = crypt.crypt(password,salt)
    def _get_password_plain(self):
        return self._password_plain
    password_plain = property(_get_password_plain, _set_password_plain)

    @classmethod
    def generatePassword(cls):
        return ''.join(random.sample(string.letters+string.digits, 12))

    @classmethod
    def getExpireDate(cls):
        return datetime.today()+ timedelta(hours=cls.ADDITIONAL_HOURS)

    def _get_is_notify(self):
        if (self.status == 'A' and self.isExpired == True  and self.confirmed == True) or (self.known_partners is not None):
            return True
        else: 
            return False
    isNotify = property(_get_is_notify)

    def _get_expire_in(self):
        if not self.expire: 
            return 'N/A'
        timediff = self.expire - datetime.today()
        return timediff

    expire_in = property(_get_expire_in)

    def _get_is_expired(self):
        if self.expire:
            if self.expire < datetime.today():
                return True
        return False
   
    isExpired = property(_get_is_expired)

    @before_insert
    def _auto_confirm(self):
        
        if self.known_partners:
            self.confirmed = True
            self.status = 'A'

    def notify(self):
        if self.isNotify == True:
            self.sendAccountMail()
        else:
            self.notify_admin()


    def notify_admin(self):
        from email.mime.text import MIMEText
        mailtext = """%(name)s of %(organisation)s has requested a temporary FTP-Account.

Please confirm the request here: %(link)s""" % dict(name=self.realname,organisation=self.organisation,link='http://arabidopsis.gmi.oeaw.ac.at:5000/management/RequestFTPAccess/confirmlist')
        msg = MIMEText(mailtext)
        me = 'arabidopsis@gmi.oeaw.ac.at'
        msg['Subject'] = 'New Temporary FTP-Account-Request'
        msg['From'] = me
        msg['To'] = self.ADMIN_MAIL
        #return True
        # Send the message via our own SMTP server, but don't include the
        # envelope header.
        s = smtplib.SMTP()
        s.connect()
        s.sendmail(me, [self.email], msg.as_string())
        s.quit()


    def sendAccountMail(self):

        from email.mime.text import MIMEText
        mailtext = """Dear Mr/Mrs %(name)s! 

You requested a Temporary sftp-Account for the Arabidopsis-Server.

You can connect to the Arabidopsis sftp-Server using following account information:
Login:       %(username)s
Password:  %(password)s
Server: %(server)s | %(direct_link)s

The account is active until %(expire)s

You can use FileZilla or any other sftp capable application to transfer data from and to the server. 

Regards
Norborg-Group""" % dict(name=self.realname,username=self.username,password=self.password_plain,server = 'arabidopsis.gmi.oeaw.ac.at:22',direct_link='sftp://%s:%s@arabidopsis.gmi.oeaw.ac.at:22' % (self.username,self.password_plain),expire = self.expire)
        msg = MIMEText(mailtext)
        me = 'arabidopsis@gmi.oeaw.ac.at'
        msg['Subject'] = 'Arabidopsis Temporary sFTP-Account'
        msg['From'] = me
        msg['To'] = self.email
        #return True
        # Send the message via our own SMTP server, but don't include the
        # envelope header.
        s = smtplib.SMTP()
        s.connect()
        s.sendmail(me, [self.email], msg.as_string())
        s.quit()


    def confirm(self):
        if not self.known_partners:
            self.known_partners = known_partners.get_by(email = self.email)
        self.confirmed = True
        if not self.known_partners: 
            self.known_partners = known_partners(email=self.email)
        self.expires = self.getExpireDate()
        self.status = 'A'
        self.password_plain = self.generatePassword()
        __session__.flush()


class known_partners(Entity):
    id = Field(Integer, primary_key=True)
    email = Field(String(100), required=True)
    created  = Field(DateTime,default=datetime.now)
    users = OneToMany('sys_users')
    using_options(tablename='known_partners')
    using_table_options(mysql_engine='InnoDB')




