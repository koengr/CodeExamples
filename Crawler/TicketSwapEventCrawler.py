from bs4 import BeautifulSoup as bs
import requests
import pandas as pd
import os
import time
import re

''' 

This is a simple crawler that checks a list of Ticketswap pages, and stores the current ticket prices for later reference. 

NB: THIS IS NOT A BOT THAT AUTOMATICALLY BUY TICKETS. 


Conventions (defined by example):
    url="https://www.ticketswap.com/event/pleinvrees-festival-2023/253d3b60-2dc6-48e1-becc-6e07be643757"
    event_id = "pleinvrees-festival-2023/253d3b60-2dc6-48e1-becc-6e07be643757"
    name = scraped from site
    date = scraped from site


DATABASES:

1. ticket_database 

Columns: [event_id, ticket_url, number_of_tickets, price ]


2. events_database

Columns: [event_id, event_url, name, date, to_monitor=True/False ]

3. crawl_database (stores all past crawls -- so that we can also observe when no tickets were available)

Columns: [ timestamp, event_url ]

'''



''' Purpose: Easy functions to open event/ticket databases, and update them. '''
class TSEventCrawler():

    ''' General information about Ticketswap HTML pages '''
    available_tickets_type = 'div'
    available_tickets_class = 'css-zmr4by emzh1320'

    sale_entry_type = 'a'

    num_tickets_type = 'h4'
    num_tickets_class = 'css-1xixvph eh8fd905'

    price_html_type = 'strong'
    price_html_class = 'css-tlq7v e1mv3lrv0'

    event_title_type = 'h1'
    event_title_class = 'css-1btrs3f ej1og8q1'

    event_date_type = 'p'
    event_date_class = 'css-mcvk2t ej1og8q8'
    ''' End general information about ticketswap HTML '''



    def __init__(self, 
            event_database_filename="event_database", 
            ticket_database_filename="ticket_database",
            crawl_database_filename="crawl_database"     ):
        
        self.event_database_filename = event_database_filename
        self.ticket_database_filename = ticket_database_filename
        self.crawl_database_filename = crawl_database_filename

        self.events = self.open_database(event_database_filename)
        self.tickets = self.open_database(ticket_database_filename)
        self.crawls = self.open_database(crawl_database_filename)


    def open_database(self, filename):
        """
        Opens a database stored in a pickle file.
        
        Parameters:
        - filename: str, the filename of the pickle file to be opened
        
        Returns:
        - A Pandas dataframe containing the data stored in the pickle file.
        
        If the file does not exist or is empty, an empty Pandas dataframe is returned.
        """

        try:
            return pd.read_pickle(filename) if os.path.isfile(filename) else pd.DataFrame()
        except EOFError:
            print(f"The database '{filename}' seems too small. Creating new empty database")
            return pd.DataFrame()


    def save(self):
        self.events.to_pickle(self.event_database_filename)
        self.tickets.to_pickle(self.ticket_database_filename)
        self.crawls.to_pickle(self.crawl_database_filename)


    ''' Open the ticketswap page through the provided URL. Return the corresponding BeautifulSoup object '''
    def open_ticketswap_url(self, url):
        headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64" }
        r = requests.get(url, headers=headers)
        return bs(r.text, 'html.parser')

    ''' Update the Tickets database with with a recording of all tickets of events that are monitored     '''
    def update(self, save=False, verbose=True):
        url_list = self.events[ self.events["To monitor"] ]["url"]
        for url in url_list: 
            if verbose: 
                print('Checking for event ' + str(url), end='   ')
            self.update_tickets_from_url( url, save=False )
            if verbose:
                print('Done!')

            # After each check, build in some safety time to avoid server bans. 
            time.sleep(0.4)

        if save :
            self.save()


    ''' 
    =============== TICKET EDITING FUNCTIONS ===================
    
    '''


    def update_tickets_from_url(self, event_url, save = False ):
        """
        Updates the ticket database with current data from a given event URL.
        
        Parameters:
        - event_url: str, the ticketswap event URL
        - save: bool, optional, default=False, indicates whether to save the updated ticket database
        
        Returns:
        - None
        
        This function fetches ticket data from the provided event URL and appends it to the current ticket database.
        It also updates the crawl database with a record of the crawl, using the earliest timestamp of the new ticket data.
        If the save parameter is set to True, the updated ticket and crawl databases are saved to pickle files.
        """
        
        # Append all currently available ticket to the ticket databse
        (new_db, timestamp) = self.ticket_database_from_url(event_url)
        self.tickets = pd.concat( [self.tickets, new_db] )

        # Register the crawl (using earliest timestamp of 'new_db')
        new_crawl = pd.DataFrame({
                'Timestamp': [timestamp],
                'Event URL': [event_url],
        })
        self.crawls = pd.concat( [self.crawls, new_crawl] )

        # Save the databases 
        if save : 
            self.save()




    def ticket_database_from_url(self, event_url):
        ''' 
        Create a Pandas Database with all the available tickets at the given URL
        '''


        '''
         Each page has 3 divs with the style of 'available tickets': the AVAILABLE ones, an EMPTY div, and the SOLD tickets. 
         The soup.find function simply returns the first div with this style. 

         Within these, we look for the head of individual sellers 

        For each seller, we want to extract three things:
        - the ticket URL
        - the number of tickets
        - the ticket price (only if it is in euros)
        '''

        # Load the page
        soup = self.open_ticketswap_url(event_url)
        timestamp = pd.Timestamp.now()

        # Search for available tickets. If none, then abort. 
        available_tickets_html = soup.find(self.available_tickets_type, {"class": self.available_tickets_class })
        sellers_html = available_tickets_html.find_all(self.sale_entry_type)

        if len(sellers_html) == None :
            return ( None, timestamp )

        # If okay, read out the ticket information.
        urls = []
        num_tickets = []
        prices = []

        for seller in sellers_html:
            price = self.get_entry_price( seller )
            if price is not None:
                urls.append( self.get_entry_url(seller) )
                num_tickets.append( self.get_entry_number_of_tickets(seller))
                prices.append(price)


        # To return: create a new Pandas dataframe that contains all the newfound information. 
        db = pd.DataFrame(
            { 
                'Timestamp': timestamp,
                'Event URL': event_url,
                'Ticket URL': urls,
                'Price': prices,
                'Num tickets': num_tickets
            }
        )
        return( (db, timestamp) )

    def get_entry_url(self, entry_html ):
        return entry_html['href']

    def get_entry_number_of_tickets(self, entry_html ):
        num_tickets_html = entry_html.find(self.num_tickets_type, {"class": self.num_tickets_class })
        num_tickets = int(num_tickets_html.contents[0].split()[0])
        return(num_tickets)

    def get_entry_price(self, entry_html ):

        # Find the string that contains the ticket price. 
        price_html = entry_html.find(self.price_html_type, {"class": self.price_html_class}).contents[0]

        # If the price is not in Euros, we don't handle the result for now. 
        if len( re.findall(r'€', price_html) ) == 0 : 
            return None
                
        # To be able to work with both decimal commas and points (and other garbage in the stirng): 
        # simply extract two numbers using regex, and assume these are (whole euros) + (cents). 
        pattern = r'\d+'
        price = re.findall( pattern, price_html )
        assert( len(price) == 2 )
        price = int(price[0]) + int(price[1])/100 
        return price 



    ''' 
    ===============  EVENT EDITING FUNCTIONS =====================
    '''

    
    ''' Given a list of events, add these to event database. Check if it isn't already present. Create a unqiue event_id. '''
    def add_events(self, url_list, to_monitor = True, save=False):

        for url in url_list:

            if (self.events.size == 0) or (url not in self.events['url'].values):
                [name, date] = self.get_event_name_and_date(url)
                self.events = self.events.append( pd.DataFrame(
                    {
                    'url': [url],
                    'Event name': [name],
                    'Event date': [date],
                    'To monitor': [to_monitor]
                    }
                ))
                print("Succesfully added event: ", url)
            
            else : 
                print("!! Event already in database: " , url )

        if save :
            self.save()

    # ''' Creates a somewhat readable filename (without illegal symbols like '/' from a TS URL 
    #     NB: we assume that 'events' is in the url, i.e. of the form https://www.ticketswap.com/event/...	
    # '''
    # def event_id_from_url(self, url):
    #     url_as_list = url.split('/')
    #     keep_from_index = url_as_list.index('event')
    #     return '_'.join( url_as_list[keep_from_index+1:])
        

    def get_event_name_and_date(self, url):
        soup = self.open_ticketswap_url(url)
        name = soup.find(self.event_title_type, {"class": self.event_title_class}).contents[0]
        date = soup.find(self.event_date_type, {"class": self.event_date_class}).contents[0]

        return [str(name), str(date)]



if __name__ == "__main__":
    ts = TSEventCrawler()
    ts.update()
